# =============================================================
# Galleri model engine (project-specific)
#
# Annual-cycle cohort simulation over a lifetime horizon (to max_age):
#   1. NATURAL HISTORY is calibrated to SEER. For each subgroup (race/
#      ethnicity x sex x Kunzmann group) the yearly onset of preclinical
#      upper GI cancer is solved so that, without screening, symptomatic
#      presentation reproduces SEER age-specific incidence exactly, and the
#      stage at symptomatic diagnosis reproduces SEER's stage distribution
#      (esophagus and stomach separately, weighted by the subgroup's site mix
#      at each age). The same calibrated natural history drives every arm.
#   2. SCREENING happens at the start of each screening year until
#      settings$screen_stop_age:
#        Kunzmann -> EGD : score once; people flagged have an EGD (one-time by default, kunz_egd_rounds);
#        Galleri         : everyone undiagnosed is tested; a +ve esophagus/
#                          stomach signal goes through the full workup tree.
#      Screen-detected cancers enter lead-time states (no cancer deaths until
#      they would have presented), so benefit comes from stage shift only.
#   3. Cancers elsewhere found by the Galleri workup are valued as the
#      difference between detection now and symptomatic presentation later,
#      using the same natural history (backward induction).
# =============================================================

# pre_m: preclinical cancer already missed by Galleri that stays undetectable by
# Galleri (scenario p_fn_persist > 0); otherwise it behaves exactly like pre.
STATES <- c("well", "pre", "pre_m", "lead_early", "lead_late", "early", "late", "dead_cancer", "dead_other")
PRE <- c("pre", "pre_m")
DEAD <- c("dead_cancer", "dead_other")
CHANNELS <- c("cost", "qaly", "daly", "ly", "clin_dx", "clin_late", "ca_death", "oc_death", "symp_perf")
ARMS <- c("Usual care", "Kunzmann -> EGD", "Galleri")
POPULATION <- "screening"
RACE_LABELS <- c(all = "All races/ethnicities", nhw = "Non-Hispanic White", nhb = "Non-Hispanic Black",
                 nhapi = "Non-Hispanic Asian/Pacific Islander",
                 nhaian = "Non-Hispanic American Indian/Alaska Native", hisp = "Hispanic")

# ---- Baseline data and calibration targets ---------------------------------
baseline_for <- function(data, row, settings) {
  B <- data$baseline_risks
  ages <- row$start_age:(settings$max_age - 1)
  d <- B[B$race_eth == row$race_eth & B$sex == row$sex, ]
  d <- d[match(ages, d$age_from), ]
  if (anyNA(d$age_from)) stop(sprintf("no baseline rows for %s %s age %d-%d", row$race_eth, row$sex,
                                     min(ages), max(ages)))
  d
}

# Kunzmann score >= 8, by sex (sex is one of the score's items, so the share
# flagged and the share of EAC captured differ by sex; RRs are within-sex).
kunz_share <- function(settings, sex) settings[[paste0("kunz_share_", sex)]]
kunz_rr <- function(risk_group, settings, sex) {
  s <- kunz_share(settings, sex); k <- settings[[paste0("kunz_capture_", sex)]]
  switch(risk_group, kunzmann_high = k / s, kunzmann_low = (1 - k) / (1 - s), 1)
}

# SEER clinical incidence of upper GI cancer for this subgroup (per person-year).
# The Kunzmann relative risk acts on the adenocarcinoma share of esophageal cancer.
target_ugi <- function(curve, row, settings) {
  rr <- kunz_rr(row$risk_group, settings, row$sex)
  curve$inc_eso * (curve$eac_share * rr + 1 - curve$eac_share) + curve$inc_stom
}

# SEER stage at symptomatic diagnosis for this subgroup, by age: share
# regional/distant, weighted by the esophagus vs stomach mix of the target.
# or_late_clin (default 1) shifts the odds for probabilistic analysis.
target_late <- function(curve, row, settings, p) {
  rr <- kunz_rr(row$risk_group, settings, row$sex)
  eso <- curve$inc_eso * (curve$eac_share * rr + 1 - curve$eac_share); stom <- curve$inc_stom
  pl <- (eso * curve$late_eso + stom * curve$late_stom) / (eso + stom)
  or <- p$or_late_clin %||% 1
  pl * or / (1 - pl + pl * or)
}

transition_matrix <- function(mo, onset, p, pl) {
  pp <- 1 - exp(-1 / p$sojourn_years); s <- 1 - mo
  n <- length(STATES)
  P <- matrix(0, n, n, dimnames = list(STATES, STATES))
  P[, "dead_other"] <- mo
  P["well", c("well", "pre")] <- s * c(1 - onset, onset)
  P["pre", c("pre", "early", "late")] <- s * c(1 - pp, pp * (1 - pl), pp * pl)
  P["pre_m", c("pre_m", "early", "late")] <- s * c(1 - pp, pp * (1 - pl), pp * pl)
  P["lead_early", c("lead_early", "early")] <- s * c(1 - pp, pp)
  P["lead_late", c("lead_late", "late")] <- s * c(1 - pp, pp)
  P["early", c("early", "dead_cancer")] <- s * c(1 - p$m_early, p$m_early)
  P["late", c("late", "dead_cancer")] <- s * c(1 - p$m_late, p$m_late)
  P[DEAD, ] <- 0; P["dead_cancer", "dead_cancer"] <- 1; P["dead_other", "dead_other"] <- 1
  P
}

# Solve yearly onset so that modelled symptomatic incidence per living person
# equals the SEER target at every age (closed form, one year at a time).
# Returns the onset vector, the starting state (preclinical prevalence implied
# by the target) and the fit.
calibrate <- function(curve, target, p, pl) {
  n <- nrow(curve); pp <- 1 - exp(-1 / p$sojourn_years); s <- 1 - curve$mort_other
  x <- setNames(numeric(length(STATES)), STATES)
  x["pre"] <- min(0.5, target[1] / (s[1] * pp)); x["well"] <- 1 - x["pre"]
  x0 <- x; onset <- numeric(n); clipped <- 0
  for (t in seq_len(n)) {
    if (t < n) {
      P0 <- transition_matrix(curve$mort_other[t], 0, p, pl[t])
      nxt <- drop(x %*% P0)                       # next state with no new onset
      alive_next <- sum(nxt[!STATES %in% DEAD])
      pre_goal <- target[t + 1] * alive_next / (s[t + 1] * pp)
      o <- (pre_goal - nxt[["pre"]]) / max(1e-12, nxt[["well"]])
      if (o < 0 || o > 1) clipped <- clipped + 1
      onset[t] <- min(1, max(0, o))
    } else onset[t] <- onset[max(1, t - 1)]
    x <- drop(x %*% transition_matrix(curve$mort_other[t], onset[t], p, pl[t]))
  }
  list(onset = onset, x0 = x0, clipped = clipped, pl = pl)
}

# ---- Rewards per cycle (as before; onset now comes from calibration) --------
cycle_rewards <- function(P, curve, t, p) {
  alive <- setNames(ceakit::alive_fraction(P, DEAD), STATES)
  live <- !STATES %in% DEAD
  R <- matrix(0, length(STATES), length(CHANNELS), dimnames = list(STATES, CHANNELS))
  u <- c(well = 1, pre = 1, pre_m = 1, lead_early = p$u_mult_early, lead_late = p$u_mult_late,
         early = p$u_mult_early, late = p$u_mult_late, dead_cancer = 0, dead_other = 0)
  dw <- c(well = 0, pre = 0, pre_m = 0, lead_early = p$dw_early, lead_late = p$dw_late,
          early = p$dw_early, late = p$dw_late, dead_cancer = 0, dead_other = 0)
  R[, "ly"] <- alive
  R[, "qaly"] <- alive * curve$utility_well[t] * u[STATES]
  R["lead_early", "cost"] <- alive[["lead_early"]] * p$c_fu_early
  R["lead_late", "cost"]  <- alive[["lead_late"]] * p$c_fu_late
  R["early", "cost"] <- alive[["early"]] * p$c_fu_early + P["early", "dead_cancer"] * p$c_eol
  R["late", "cost"]  <- alive[["late"]] * p$c_fu_late  + P["late", "dead_cancer"]  * p$c_eol
  R[, "daly"] <- alive * dw[STATES] + P[, "dead_cancer"] * curve$daly_ref_le[t] * live
  # symptomatic diagnosis: workup (incl. an EGD that can perforate) and first-year treatment
  perf <- p$p_perf_egd %||% 0; harm_perf <- perf * (1 - (p$p_death_perf %||% 0)) * (p$d_sae %||% 0)
  for (st in PRE) {
    dx <- P[st, "early"] + P[st, "late"]
    R[st, "cost"] <- P[st, "early"] * (p$c_symp_dx + p$c_tx_first_early) +
                     P[st, "late"]  * (p$c_symp_dx + p$c_tx_first_late) + dx * perf * (p$c_perf %||% 0)
    R[st, "qaly"] <- R[st, "qaly"] - dx * harm_perf
    R[st, "daly"] <- R[st, "daly"] + dx * harm_perf
    R[st, "clin_dx"] <- dx
    R[st, "clin_late"] <- P[st, "late"]
    R[st, "symp_perf"] <- dx * perf
  }
  R[, "ca_death"] <- P[, "dead_cancer"] * live
  R[, "oc_death"] <- P[, "dead_other"] * live
  R
}

disc_vec <- function(settings, t) {   # discount factors for cycle t (t = 1 is year 0)
  r <- c(cost = settings$disc_cost, qaly = settings$disc_qaly, daly = settings$disc_daly)
  d <- setNames(rep(1, length(CHANNELS)), CHANNELS)
  d[names(r)] <- (1 + r)^-(t - 1)
  d
}

# Backward induction: value (discounted to the start of cycle t) of being in
# each state at the start of cycle t, under usual care. Used to value cancers
# elsewhere that the Galleri workup finds.
values_by_cycle <- function(Ps, Rs, settings) {
  n <- length(Ps); V <- vector("list", n + 1)
  V[[n + 1]] <- matrix(0, length(STATES), length(CHANNELS), dimnames = list(STATES, CHANNELS))
  d1 <- disc_vec(settings, 2)
  for (t in n:1) V[[t]] <- Rs[[t]] + (Ps[[t]] %*% V[[t + 1]]) * rep(d1, each = length(STATES))
  V[seq_len(n)]
}

# ---- Galleri +ve workup (revised Nov 2026) ----------------------------------
# Step 1, the EGD, is applied in galleri_round(): it is skipped inside the hold after
# a negative EGD (egd_hold_years). Step 2, the rest of the workup after a negative
# or skipped EGD, is the tree below. It depends on pancreatobiliary history:
# 1 = none, 2 = prior negative EUS, 3 = prior negative laparoscopy.
PB <- c("none", "eus", "lap")

rest_tree <- function(p, pb, edges_only = FALSE) {
  M <- p$p_pb_malig; se <- p$eus_tissue_sens; ps <- p$p_eus_samp
  E <- list(); L <- list()
  edge <- function(f, t, pr) E[[length(E) + 1]] <<- data.frame(from = f, to = t, prob = pr)
  leaf <- function(id, found, missed, flag) L[[length(L) + 1]] <<- data.frame(leaf = id, found = found, missed = missed, flag = flag)
  edge("start", "ent", p$p_ent_route); edge("start", "pet", 1 - p$p_ent_route)
  edge("ent", "ent_pos", p$p_ent_pos); edge("ent", "ent_neg", 1 - p$p_ent_pos)
  leaf("ent_pos", 1, 0, pb); leaf("ent_neg", 0, p$p_miss, pb)
  edge("pet", "pb", p$p_pet_pb); edge("pet", "incid", p$p_pet_incid); edge("pet", "pet_neg", 1 - p$p_pet_pb - p$p_pet_incid)
  leaf("incid", p$p_incid_cancer, 0, pb); leaf("pet_neg", 0, p$p_miss, pb)
  lap <- function(from, pr, m, tag) {                     # diagnostic laparoscopy (finds malignancy present)
    edge(from, paste0("lap.", tag), pr)
    edge(paste0("lap.", tag), paste0("lappos.", tag), m); edge(paste0("lap.", tag), paste0("lapneg.", tag), 1 - m)
    leaf(paste0("lappos.", tag), 1, 0, 3); leaf(paste0("lapneg.", tag), 0, 0, 3)
  }
  ww <- function(from, pr, m, tag, flag) {                # watchful waiting: cancers present are missed for now
    edge(from, paste0("ww.", tag), pr); leaf(paste0("ww.", tag), 0, m, flag)
  }
  eus <- function(from, pr, m, tag) {                     # EUS; FNA/FNB a coin toss
    n <- paste0("eus.", tag); edge(from, n, pr)
    edge(n, paste0("samp.", tag), ps); edge(n, paste0("nosamp.", tag), 1 - ps)
    edge(paste0("samp.", tag), paste0("euspos.", tag), m * se); edge(paste0("samp.", tag), paste0("eusneg.", tag), 1 - m * se)
    leaf(paste0("euspos.", tag), 1, 0, 2)
    list(neg = paste0("eusneg.", tag), m_neg = m * (1 - se) / (1 - m * se), nos = paste0("nosamp.", tag), m_nos = m)
  }
  if (pb == 1) {                                          # first finding: EUS, then tumor board
    r <- eus("pb", 1, M, "a")
    for (z in list(list(r$neg, r$m_neg, "n"), list(r$nos, r$m_nos, "s"))) {
      tb <- paste0("tb.", z[[3]]); edge(z[[1]], tb, 1)
      lap(tb, p$p_tb_lap, z[[2]], z[[3]]); ww(tb, 1 - p$p_tb_lap, z[[2]], z[[3]], 2)
    }
  } else {                                                # repeat finding: laparoscopy less likely
    pl <- p$p_tb_lap * if (pb == 2) p$r_lap_after_eus else p$r_lap_after_lap
    pe <- if (pb == 2) p$p_reeus_after_eus else p$p_reeus_after_lap
    lap("pb", pl, M, "r")
    r <- eus("pb", pe, M, "r")
    ww(r$neg, 1, r$m_neg, "rn", pb); ww(r$nos, 1, r$m_nos, "rs", pb)
    ww("pb", 1 - pl - pe, M, "r0", pb)
  }
  if (isTRUE(edges_only)) return(do.call(rbind, E))
  d <- ceakit::tree_paths(do.call(rbind, E))
  d <- merge(d, do.call(rbind, L), by = "leaf", sort = FALSE)
  tok <- strsplit(d$path, " > ", fixed = TRUE)
  cnt <- function(type) vapply(tok, function(z) sum(sub("\\..*$", "", z) == type), 0)
  for (k in c("ent", "pet", "eus", "samp", "tb", "lap", "ww", "incid")) d[[paste0("n_", k)]] <- cnt(k)
  d$free <- pmax(0, 1 - d$found - d$missed)
  d$cost <- d$n_ent * p$c_ent + d$n_pet * p$c_pet + d$n_eus * p$c_eus + d$n_samp * p$c_eus_samp + d$n_tb * p$c_tb +
            d$n_lap * p$c_lap + d$n_ww * p$c_surv + d$n_incid * (p$c_incid_fu %||% 0)
  d$days <- d$n_ent * p$t_ent + d$n_pet * p$t_pet + d$n_eus * p$t_eus + d$n_tb * p$t_tb + d$n_lap * p$t_lap
  d$sae <- d$n_ent * p$sae_ent + d$n_eus * p$sae_eus + d$n_samp * p$sae_eus_samp + d$n_lap * p$p_lap_cx
  d$rad <- d$n_pet * p$msv_pet / 1000 * p$risk_per_sv
  d$invasive <- d$n_eus + d$n_samp + d$n_lap
  d$dfp <- ifelse(startsWith(d$leaf, "ent"), p$dfp_endo, ifelse(d$leaf == "pet_neg", p$dfp_pet, p$dfp_ext))
  d
}

# Averages over the rest of the workup: for people with another cancer (all paths,
# weighted by probability) and for people with no cancer outside the upper GI tract
# (paths weighted by their cancer-free share), plus where the latter's
# pancreatobiliary history ends up.
REST_COLS <- c("n_ent", "n_pet", "n_eus", "n_samp", "n_tb", "n_lap", "n_ww", "n_incid", "cost", "days", "sae", "rad", "invasive", "dfp")
rest_summary <- function(p, pb) {
  d <- rest_tree(p, pb)
  a <- colSums(d[REST_COLS] * d$prob)
  w <- d$prob * d$free; f <- colSums(d[REST_COLS] * w) / sum(w)
  to <- vapply(1:3, function(b) sum(w[d$flag == b]) / sum(w), 0)
  list(all = a, free = f, det_o = sum(d$prob * d$found) / sum(d$prob * (d$found + d$missed)), to = to, paths = d)
}
workup_summary <- function(p) list(rest = lapply(1:3, function(b) rest_summary(p, b)))
rest_tree_edges <- function(p, pb = 1) rest_tree(p, pb, edges_only = TRUE)

# Probability that a preclinical cancer present at the start of cycle t would
# NEVER be diagnosed with symptoms within the horizon (other-cause death first, or
# still preclinical at max_age). A screen-detected cancer with this outcome is an
# overdiagnosis. Lead-time states progress and die exactly like 'pre', so this is
# exact for upper GI cancers; other-site cancers found share the same natural history.
overdx_prob <- function(curve, cal, p) {
  n <- nrow(curve); q <- numeric(n + 1); q[n + 1] <- 1
  for (t in n:1) {
    P <- transition_matrix(curve$mort_other[t], cal$onset[t], p, cal$pl[t])
    q[t] <- P["pre", "dead_other"] + P["pre", "pre"] * q[t + 1]
  }
  q[seq_len(n)]
}

# ---- One cohort row, one arm: annual simulation -----------------------------
EVENTS <- c("kunz", "tests", "investigated", "egd", "ent", "pet", "eus", "samp", "tb", "lap", "surv", "invasive",
            "invasive_free", "free_inv", "sae", "sae_free", "proc_deaths", "rad_ca", "days", "found_ugi",
            "found_other", "missed_inv", "late_screen", "rounds",
            "tp", "fp", "fp_other_ca", "tn", "fn", "fn_symp", "egd_pos", "egd_neg", "perf", "perf_deaths",
            "overdx_ugi", "overdx_other", "egd_free", "ent_free", "pet_free", "eus_free", "samp_free", "lap_free", "surv_free",
            "perf_free", "deaths_free", "rad_free", "days_free", "incid", "incid_free", "fn_persist",
            "worry", "sae_other_free", "egd_skipped", "found_pet_ugi", "hidden_new", "missed_visible", "repeat_fp")

# History tags carried by people not yet diagnosed in the Galleri arm (Nov 2026):
#   k  = earlier esophagus/stomach false positives (0, 1, 2+): each person's false-positive
#        propensity is beta(m kappa, (1 - m) kappa), updated by Bayes with every result;
#   h  = annual rounds left in the EGD hold after a negative EGD (0 = EGD if positive);
#   b  = pancreatobiliary history (none, prior negative EUS, prior negative laparoscopy).
# Health columns: well, pre (preclinical, visible on EGD), iv (preclinical, not yet
# visible), pm (missed by Galleri and staying undetectable: scenario p_fn_persist).
HCOLS <- c("well", "pre", "iv", "pm")
tag_grid <- function(H) expand.grid(k = 0:2, h = 0:H, b = 1:3)
tag_id <- function(T, k, h, b) which(T$k == k & T$h == h & T$b == b)

# One screening round of the Galleri arm, all history tags at once (vectorised).
tag_idx <- function(k, h, b, H) 1 + k + 3 * h + 3 * (H + 1) * (b - 1)   # row of tag_grid(H)
acc_rows <- function(z, j, v) { r <- rowsum(v, j); i <- as.integer(rownames(r)); z[i] <- z[i] + r[, 1]; z }
galleri_round <- function(X, T, t, nprev, curve, p, settings, wk) {
  H <- max(T$h); interval <- max(1, settings$screen_interval); n <- nrow(T)
  rep_m <- if (nprev > 0) p$sens_g_rep_ratio %||% 1 else 1
  m <- p$fpr_g_es; kap <- p$fp_kappa %||% 1e6
  sE <- 1 - p$p_perf_egd * p$p_death_perf; fi <- p$f_invisible %||% 0; se_egd <- p$sens_egd
  reach_pet <- 1 - p$p_ent_route; ppu <- p$p_pet_ugi %||% 0
  Fm <- do.call(rbind, lapply(wk$rest, `[[`, "free"))[T$b, , drop = FALSE]   # per tag: cancer-free averages
  Am <- do.call(rbind, lapply(wk$rest, `[[`, "all"))[T$b, , drop = FALSE]    # per tag: other-cancer averages
  det <- vapply(wk$rest, `[[`, 0, "det_o")[T$b]
  TO <- do.call(rbind, lapply(wk$rest, `[[`, "to"))[T$b, , drop = FALSE]
  W <- X[, "well"]; U <- X[, "pre"]; I <- X[, "iv"]; Mm <- X[, "pm"]
  k <- T$k; h <- T$h; b <- T$b
  theta <- (m * kap + k) / (kap + nprev)                           # posterior mean false-positive propensity
  pool_o <- W * curve$inc_other[t] * (if (t == 1) p$sojourn_years else interval)
  a_u <- U * p$sens_g_ugi * rep_m; a_i <- I * p$sens_g_ugi * rep_m
  a_o <- pool_o * p$sens_g_other * rep_m; a_f <- (W - pool_o) * theta
  gp <- a_u + a_i + a_o + a_f
  egd <- h == 0; s_e <- ifelse(egd, sE, 1)
  fu1 <- ifelse(egd, a_u * se_egd, 0)
  mv <- ifelse(egd, a_u * (1 - se_egd) * (1 - fi), a_u)          # visible, missed (or EGD skipped)
  mi <- ifelse(egd, a_u * (1 - se_egd) * fi, 0) + a_i              # not yet visible
  fu2 <- mv * s_e * reach_pet * ppu                                # seen on PET -> targeted EGD
  fo <- a_o * s_e * det
  n_egd <- ifelse(egd, gp, 0) + fu2
  rest_free <- (a_f + mv + mi) * s_e; rest_o <- a_o * s_e
  R <- function(col) sum(rest_free * Fm[, col] + rest_o * Am[, col])
  perf <- n_egd * p$p_perf_egd; sae_o_v <- rest_free * Fm[, "sae"] + rest_o * Am[, "sae"]
  pd <- perf * p$p_death_perf + sae_o_v * p$p_death_sae
  ff <- a_f * s_e; pf <- ifelse(egd, a_f * p$p_perf_egd, 0); af_e <- ifelse(egd, a_f, 0)
  FF <- function(col) sum(ff * Fm[, col])
  e <- setNames(as.list(numeric(length(EVENTS))), EVENTS)
  e$tests <- sum(W + U + I + Mm); e$investigated <- sum(gp); e$free_inv <- sum(a_f)
  e$found_ugi <- sum(fu1 + fu2); e$found_other <- sum(fo); e$found_pet_ugi <- sum(fu2)
  e$missed_inv <- sum(a_u + a_i + a_o - fu1 - fu2 - fo); e$egd_skipped <- sum(gp[!egd])
  e$egd <- sum(n_egd); e$egd_pos <- sum(fu1); e$egd_neg <- sum(gp[egd]) - sum(fu1)
  for (kk in c("ent", "pet", "eus", "samp", "tb", "lap")) e[[kk]] <- R(paste0("n_", kk))
  e$surv <- R("n_ww"); e$incid <- R("n_incid")
  e$invasive <- sum(n_egd) + R("invasive"); e$perf <- sum(perf); e$perf_deaths <- sum(perf) * p$p_death_perf
  e$sae <- sum(perf) + sum(sae_o_v); e$rad_ca <- R("rad")
  e$days <- sum(gp[egd]) * p$t_egd + R("days")
  e$tp <- sum(a_u + a_i); e$fp <- sum(a_o + a_f); e$fp_other_ca <- sum(a_o); e$tn <- sum(W - a_o - a_f)
  e$fn <- sum(U + I - a_u - a_i + Mm)
  e$hidden_new <- sum(mi); e$missed_visible <- sum(mv - fu2); e$repeat_fp <- sum(a_f[k > 0])
  e$egd_free <- sum(af_e); e$perf_free <- sum(pf)
  for (kk in c("ent", "pet", "eus", "samp", "incid", "lap")) e[[paste0(kk, "_free")]] <- FF(paste0("n_", kk))
  e$surv_free <- FF("n_ww")
  e$invasive_free <- sum(af_e) + FF("invasive")
  e$sae_other_free <- FF("sae"); e$sae_free <- sum(pf) + FF("sae")
  e$deaths_free <- sum(pf) * p$p_death_perf + FF("sae") * p$p_death_sae
  e$rad_free <- FF("rad"); e$days_free <- sum(af_e) * p$t_egd + FF("days")
  e$worry <- sum(a_f) * (p$p_worry %||% 0)
  persist <- (U - a_u) * (p$p_fn_persist %||% 0); e$fn_persist <- sum(persist)
  q_path <- FF("dfp")
  fp_cost <- sum(af_e) * (p$c_egd + p$p_perf_egd * p$c_perf) + FF("cost") + FF("sae") * p$c_sae
  cost_t <- sum(W + U + I + Mm) * p$c_galleri * (p$galleri_cost_share %||% 1) + sum(n_egd) * p$c_egd + R("cost")
  # where everyone goes next round
  h_after <- pmax(h - 1, 0); h_pos <- ifelse(egd, H, h_after)
  stay <- tag_idx(k, h_after, b, H)
  Xn <- X * 0
  Xn[, "well"] <- acc_rows(Xn[, "well"], stay, W - a_o - a_f - pd)
  Xn[, "pre"]  <- acc_rows(Xn[, "pre"], stay, U - a_u - persist)
  Xn[, "iv"]   <- acc_rows(Xn[, "iv"], stay, I - a_i)
  Xn[, "pm"]   <- acc_rows(Xn[, "pm"], stay, Mm + persist)
  for (b2 in 1:3) {                                                # positives not diagnosed: hold and history updated
    w2 <- TO[, b2]
    j_fp <- tag_idx(pmin(k + 1, 2), h_pos, b2, H); j <- tag_idx(k, h_pos, b2, H)
    Xn[, "well"] <- acc_rows(Xn[, "well"], j_fp, a_f * w2)
    Xn[, "well"] <- acc_rows(Xn[, "well"], j, a_o * w2)
    Xn[, "pre"]  <- acc_rows(Xn[, "pre"], j, (mv - fu2) * w2)
    Xn[, "iv"]   <- acc_rows(Xn[, "iv"], j, mi * w2)
  }
  list(X = Xn, e = lapply(e, function(v) unname(v)), q_path = unname(q_path), fp_cost = unname(fp_cost),
       cost_t = unname(cost_t), sae_other = unname(sum(sae_o_v)), pd = unname(sum(pd)))
}

# Natural history for the tagged people (the same transitions as transition_matrix);
# returns the next tagged matrix and the flows into diagnosed and dead states.
tagged_step <- function(X, T, P, onset, p) {
  s <- 1 - P["well", "dead_other"]; pp <- P["pre", "early"] + P["pre", "late"]
  pl <- if (pp > 0) P["pre", "late"] / pp else 0
  v <- 1 - exp(-1 / (p$t_visible %||% 1)); hr <- ifelse(T$k > 0, p$hr_ca_after_fp %||% 1, 1)
  sp <- if (s > 0) pp / s else 0                                 # presentation among survivors
  on <- pmin(1, onset * hr)
  clin <- (X[, "pre"] + X[, "iv"] + X[, "pm"]) * pp
  Xn <- X * 0
  Xn[, "well"] <- X[, "well"] * s * (1 - on)
  Xn[, "pre"] <- X[, "well"] * s * on + X[, "pre"] * (s - pp) + X[, "iv"] * (s - pp) * v
  Xn[, "iv"] <- X[, "iv"] * (s - pp) * (1 - v)
  Xn[, "pm"] <- X[, "pm"] * (s - pp)
  list(X = Xn, early = unname(sum(clin) * (1 - pl)), late = unname(sum(clin) * pl), dead_other = unname(sum(X) * (1 - s)))
}

agg_state <- function(X, y) { x <- y; x[["well"]] <- sum(X[, "well"]); x[["pre"]] <- sum(X[, "pre"]) + sum(X[, "iv"])
  x[["pre_m"]] <- sum(X[, "pm"]); x }

simulate_arm <- function(arm, row, curve, cal, p, settings, wk, V, q_od = overdx_prob(curve, cal, p)) {
  n <- nrow(curve); pls <- p$p_late_screen; x <- cal$x0
  cost_fp <- 0; qaly_fp <- 0
  fpq <- c(pathway = 0, worry = 0, complications = 0, deaths = 0, radiation = 0)   # QALYs lost by cancer-free people worked up
  fpd <- fpq; names(fpd)[1] <- "anxiety"                                            # DALYs, same people
  ev <- setNames(numeric(length(EVENTS)), EVENTS)
  acc <- setNames(numeric(length(CHANNELS)), CHANNELS)
  cost0 <- 0; q_harm <- 0; d_harm <- 0
  trace <- matrix(0, n + 1, length(STATES), dimnames = list(NULL, STATES)); trace[1, ] <- x
  ev_t <- matrix(0, n, length(EVENTS), dimnames = list(NULL, EVENTS))        # screening events by cycle
  acc_t <- matrix(0, n, length(CHANNELS), dimnames = list(NULL, CHANNELS))   # undiscounted channel flows by cycle
  interval <- max(1, settings$screen_interval)
  is_g <- arm == "Galleri"
  if (is_g) {                                                   # tagged people (history), all start untagged
    T <- tag_grid(max(0, round(p$egd_hold_years %||% 0)))
    X <- matrix(0, nrow(T), length(HCOLS), dimnames = list(NULL, HCOLS))
    X[tag_idx(0, 0, 1, max(T$h)), c("well", "pre", "pm")] <- c(x[["well"]], x[["pre"]], x[["pre_m"]])
  }
  for (t in seq_len(n)) {
    age <- row$start_age + t - 1
    df <- disc_vec(settings, t)
    P <- transition_matrix(curve$mort_other[t], cal$onset[t], p, cal$pl[t])
    if (arm == "Kunzmann -> EGD" && t == 1) { ev["kunz"] <- 1; cost0 <- cost0 + df[["cost"]] * p$c_kunz }
    screen <- arm != "Usual care" && age <= settings$screen_stop_age && (t - 1) %% interval == 0 &&
              (is_g || (isTRUE(row$kunzmann_scoped == 1) &&
                        ev[["rounds"]] < (settings$kunz_egd_rounds %||% Inf)))   # Kunzmann: one-time EGD by default
    if (screen) {
      if (is_g) {
        W <- sum(X[, "well"])
        g <- galleri_round(X, T, t, ev[["rounds"]], curve, p, settings, wk)
        X <- g$X; e <- g$e; sae_other <- g$sae_other; fp_cost <- g$fp_cost; cost_t <- g$cost_t
        pd <- g$pd
      } else {
        W <- x[["well"]]; U <- x[["pre"]]; e <- setNames(as.list(numeric(length(EVENTS))), EVENTS)
        egd <- W + U
        e$found_ugi <- U * p$sens_egd; e$missed_inv <- U - e$found_ugi
        e$investigated <- egd; e$free_inv <- W
        e$egd <- egd; e$invasive <- egd; e$invasive_free <- W
        e$perf <- egd * p$p_perf_egd; e$perf_deaths <- e$perf * p$p_death_perf; sae_other <- 0
        e$egd_free <- W; e$perf_free <- W * p$p_perf_egd; e$deaths_free <- e$perf_free * p$p_death_perf
        e$days_free <- W * p$t_egd
        fp_cost <- W * (p$c_egd + p$p_perf_egd * p$c_perf)
        e$sae <- e$perf; e$sae_free <- W * p$p_perf_egd; e$days <- egd * p$t_egd
        e$egd_pos <- e$found_ugi; e$egd_neg <- egd - e$found_ugi   # EGD is the screening test here
        e$tp <- e$found_ugi; e$fn <- U - e$found_ugi; e$tn <- W
        cost_t <- egd * p$c_egd
        pd <- min(W * 0.5, e$perf_deaths)
      }
      e$fn_symp <- e$fn * (P["pre", "early"] + P["pre", "late"])   # false negatives presenting with symptoms within a year
      fu <- e$found_ugi; fo <- e$found_other
      e$overdx_ugi <- fu * q_od[t]; e$overdx_other <- fo * q_od[t]
      # Burden in cancer-free people worked up. Galleri false positives: the whole workup
      # pathway (procedures and uncertainty) is valued with Matza 2024 TTO disutilities
      # (fp_tto = 1) or per procedure (d_proc); persisting worry, complications, deaths and
      # radiation cancers are added. DALYs: GBD anxiety weight x days of uncertainty, worry,
      # complications, YLL of deaths, radiation cancers. Kunzmann: per screening EGD.
      n_proc_free <- e$egd_free + e$ent_free + e$pet_free + e$eus_free + e$incid_free + e$lap_free + e$surv_free
      dw_a <- p$dw_anx %||% 0
      q_path <- if (is_g && (p$fp_tto %||% 0) == 1) g$q_path else n_proc_free * p$d_proc
      worry_w <- e$worry * (p$dur_worry %||% 0) * dw_a
      q_comp <- (e$sae_free - e$deaths_free) * p$d_sae; q_rad <- e$rad_free * p$loss_rad_ca
      fq <- c(q_path, worry_w, q_comp, e$deaths_free * V[[t]]["well", "qaly"], q_rad)
      fd <- c(if (is_g) e$days_free / 365.25 * dw_a else n_proc_free * p$d_proc, worry_w, q_comp,
              e$deaths_free * curve$daly_ref_le[t], q_rad)
      fpq <- fpq + df[["qaly"]] * fq; fpd <- fpd + df[["daly"]] * fd
      qaly_fp <- qaly_fp + df[["qaly"]] * sum(fq)
      cost_fp <- cost_fp + df[["cost"]] * fp_cost
      e$proc_deaths <- pd; e$late_screen <- (fu + fo) * pls; e$rounds <- 1
      for (k in names(e)) { ev[k] <- ev[k] + e[[k]]; ev_t[t, k] <- e[[k]] }
      # screen-detected upper GI cancers into lead time; procedure deaths out of well
      if (!is_g) { x["pre"] <- x["pre"] - fu; x["well"] <- x["well"] - pd }
      x["lead_early"] <- x["lead_early"] + fu * (1 - pls); x["lead_late"] <- x["lead_late"] + fu * pls
      x["dead_other"] <- x["dead_other"] + pd
      if (is_g) x <- agg_state(X, x)
      tx_first <- (fu + fo) * ((1 - pls) * p$c_tx_first_early + pls * p$c_tx_first_late)
      cost0 <- cost0 + df[["cost"]] * (cost_t + tx_first + e$perf * p$c_perf + sae_other * p$c_sae)
      n_proc <- e$egd + e$ent + e$pet + e$eus + e$incid + e$lap + e$surv
      # people with cancer keep the per-procedure loss; cancer-free Galleri positives use fq/fd above
      n_nf <- n_proc - (if (is_g) n_proc_free else 0)
      oth <- (e$sae - pd) * p$d_sae + e$rad_ca * p$loss_rad_ca
      qh <- n_nf * p$d_proc + oth + (if (is_g) q_path + worry_w else 0)
      dh <- n_nf * p$d_proc + oth + (if (is_g) fd[1] + worry_w else 0)
      q_harm <- q_harm + df[["qaly"]] * qh
      d_harm <- d_harm + df[["daly"]] * (dh + pd * curve$daly_ref_le[t])
      # cancers elsewhere found now vs presenting with symptoms later
      if (fo > 0) {
        v_det <- (1 - pls) * V[[t]]["lead_early", ] + pls * V[[t]]["lead_late", ]
        delta <- (v_det - V[[t]]["pre", ]) * df
        acc <- acc + fo * delta
      }
    }
    R <- cycle_rewards(P, curve, t, p)
    acc_t[t, ] <- drop(x %*% R)
    acc <- acc + acc_t[t, ] * df
    if (is_g) {
      st <- tagged_step(X, T, P, cal$onset[t], p)
      y <- x; y[c("well", "pre", "pre_m")] <- 0
      y <- drop(y %*% P)
      y[["early"]] <- y[["early"]] + st$early; y[["late"]] <- y[["late"]] + st$late
      y[["dead_other"]] <- y[["dead_other"]] + st$dead_other
      X <- st$X; x <- agg_state(X, y)
    } else x <- drop(x %*% P)
    trace[t + 1, ] <- x
  }
  list(ev = ev, acc = acc, cost0 = cost0, q_harm = q_harm, d_harm = d_harm, trace = trace,
       pre0 = cal$x0[["pre"]], cost_fp = cost_fp, qaly_fp = qaly_fp, fpq = fpq, fpd = fpd, ev_t = ev_t, acc_t = acc_t)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Outcomes per cohort row ------------------------------------------------
ADDITIVE <- c("total_cancers", "cancer_free", "kunz_scores", "galleri_tests", "screen_rounds",
  "screen_tp", "screen_fp", "screen_fp_other_cancer", "screen_tn", "screen_fn", "fn_symptomatic_1y",
  "egd_positive", "egd_negative_workup", "perforations", "perforation_deaths",
  "overdiagnosed", "overdiagnosed_ugi", "overdiagnosed_other", "egd_free", "ent_free", "pet_ct_free", "eus_free",
  "eus_sampling_free", "perforations_free", "proc_deaths_free", "dose_sv_free", "radiation_cancers_free",
  "days_free", "cost_fp", "qaly_fp", "incidental_followups", "incidental_followups_free",
  "perforations_symptomatic", "fn_persistent", "fp_worry_cases", "fp_sae_other",
  "fp_qaly_pathway", "fp_qaly_worry", "fp_qaly_complications", "fp_qaly_deaths", "fp_qaly_radiation",
  "fp_daly_anxiety", "fp_daly_worry", "fp_daly_complications", "fp_daly_deaths", "fp_daly_radiation", "daly_fp",
  "investigated", "laparoscopy", "surveillance_imaging", "laparoscopy_free", "surveillance_free",
  "egd_skipped_hold", "ugi_found_by_pet", "hidden_after_egd", "missed_visible_after_egd", "repeat_false_positives",
  "egd", "ent", "pet_ct", "eus", "eus_sampling", "tumor_board", "invasive_total", "days_x_inv",
  "cost_tests", "cost_total", "screen_detected", "found_upper_gi", "found_other",
  "late_lifetime", "clin_dx_lifetime", "ca_deaths", "oc_deaths", "ly_total", "free_investigated",
  "invasive_free", "sae_total", "sae_free", "proc_deaths", "dose_sv", "radiation_cancers",
  "missed_investigated", "qaly_total", "qaly_harms", "daly_total", "daly_harms")

# Site-specific sojourn scenario (site_sojourn = 1): the subgroup's mix of esophageal
# adenocarcinoma (sojourn_eac) and other upper GI cancer (sojourn_years), by SEER incidence.
row_params <- function(p, curve, row, settings) {
  if (!isTRUE((p$site_sojourn %||% 0) == 1)) return(p)
  rr <- kunz_rr(row$risk_group, settings, row$sex)
  eac <- curve$inc_eso * curve$eac_share * rr; tot <- target_ugi(curve, row, settings)
  w <- sum(eac) / sum(tot)
  p$sojourn_years <- w * p$sojourn_eac + (1 - w) * p$sojourn_years
  p
}

# Kunzmann sex-specific scenario (kunz_sexspec = 1): Kunzmann 2018 shares by sex
KUNZ_SEX <- list(kunz_share_male = 0.597, kunz_share_female = 0.035, kunz_capture_male = 0.881, kunz_capture_female = 0.153)
kunz_settings <- function(p, settings) {
  if (isTRUE((p$kunz_sexspec %||% 0) == 1)) settings[names(KUNZ_SEX)] <- KUNZ_SEX
  settings
}
kunz_cohort <- function(cohort, settings) {
  for (i in seq_len(nrow(cohort))) {
    k <- settings[[paste0("kunz_share_", cohort$sex[i])]]
    if (cohort$risk_group[i] == "kunzmann_high") cohort$share[i] <- k
    if (cohort$risk_group[i] == "kunzmann_low") cohort$share[i] <- 1 - k
  }
  cohort
}

prepare_row <- function(row, data, p, settings) {
  curve <- baseline_for(data, row, settings)
  p <- row_params(p, curve, row, settings)
  cal <- calibrate(curve, target_ugi(curve, row, settings), p, target_late(curve, row, settings, p))
  Ps <- lapply(seq_len(nrow(curve)), function(t) transition_matrix(curve$mort_other[t], cal$onset[t], p, cal$pl[t]))
  Rs <- lapply(seq_len(nrow(curve)), function(t) cycle_rewards(Ps[[t]], curve, t, p))
  q <- numeric(length(Ps) + 1); q[length(q)] <- 1                    # overdiagnosis probability (see overdx_prob)
  for (t in rev(seq_along(Ps))) q[t] <- Ps[[t]]["pre", "dead_other"] + Ps[[t]]["pre", "pre"] * q[t + 1]
  list(curve = curve, cal = cal, V = values_by_cycle(Ps, Rs, settings), q = q[seq_along(Ps)], p = p)
}

row_outcomes <- function(row, data, p, settings, arms, prep = prepare_row(row, data, p, settings),
                         wk = workup_summary(p)) {
  p <- prep$p
  one <- function(a) {
    s <- simulate_arm(a, row, prep$curve, prep$cal, p, settings, wk, prep$V, prep$q); e <- s$ev; l <- s$acc
    c(total_cancers = s$pre0, cancer_free = 1 - s$pre0, kunz_scores = e[["kunz"]], galleri_tests = e[["tests"]],
      screen_rounds = e[["rounds"]], screen_tp = e[["tp"]], screen_fp = e[["fp"]],
      screen_fp_other_cancer = e[["fp_other_ca"]], screen_tn = e[["tn"]], screen_fn = e[["fn"]],
      fn_symptomatic_1y = e[["fn_symp"]], egd_positive = e[["egd_pos"]], egd_negative_workup = e[["egd_neg"]],
      perforations = e[["perf"]], perforation_deaths = e[["perf_deaths"]],
      overdiagnosed = e[["overdx_ugi"]] + e[["overdx_other"]], overdiagnosed_ugi = e[["overdx_ugi"]],
      overdiagnosed_other = e[["overdx_other"]], egd_free = e[["egd_free"]], ent_free = e[["ent_free"]],
      pet_ct_free = e[["pet_free"]], eus_free = e[["eus_free"]], eus_sampling_free = e[["samp_free"]],
      perforations_free = e[["perf_free"]], proc_deaths_free = e[["deaths_free"]],
      dose_sv_free = e[["pet_free"]] * p$msv_pet / 1000, radiation_cancers_free = e[["rad_free"]],
      days_free = e[["days_free"]], cost_fp = s$cost_fp, qaly_fp = s$qaly_fp,
      incidental_followups = e[["incid"]], incidental_followups_free = e[["incid_free"]],
      perforations_symptomatic = l[["symp_perf"]], fn_persistent = e[["fn_persist"]],
      fp_worry_cases = e[["worry"]], fp_sae_other = e[["sae_other_free"]],
      fp_qaly_pathway = s$fpq[["pathway"]], fp_qaly_worry = s$fpq[["worry"]], fp_qaly_complications = s$fpq[["complications"]],
      fp_qaly_deaths = s$fpq[["deaths"]], fp_qaly_radiation = s$fpq[["radiation"]],
      fp_daly_anxiety = s$fpd[["anxiety"]], fp_daly_worry = s$fpd[["worry"]], fp_daly_complications = s$fpd[["complications"]],
      fp_daly_deaths = s$fpd[["deaths"]], fp_daly_radiation = s$fpd[["radiation"]], daly_fp = sum(s$fpd), investigated = e[["investigated"]],
      laparoscopy = e[["lap"]], surveillance_imaging = e[["surv"]], laparoscopy_free = e[["lap_free"]],
      surveillance_free = e[["surv_free"]], egd_skipped_hold = e[["egd_skipped"]], ugi_found_by_pet = e[["found_pet_ugi"]],
      hidden_after_egd = e[["hidden_new"]], missed_visible_after_egd = e[["missed_visible"]],
      repeat_false_positives = e[["repeat_fp"]], egd = e[["egd"]], ent = e[["ent"]],
      pet_ct = e[["pet"]], eus = e[["eus"]], eus_sampling = e[["samp"]], tumor_board = e[["tb"]],
      invasive_total = e[["invasive"]], days_x_inv = e[["days"]],
      cost_tests = s$cost0, cost_total = s$cost0 + l[["cost"]],
      screen_detected = e[["found_ugi"]] + e[["found_other"]], found_upper_gi = e[["found_ugi"]],
      found_other = e[["found_other"]], late_lifetime = e[["late_screen"]] + l[["clin_late"]],
      clin_dx_lifetime = l[["clin_dx"]], ca_deaths = l[["ca_death"]], oc_deaths = l[["oc_death"]],
      ly_total = l[["ly"]], free_investigated = e[["free_inv"]], invasive_free = e[["invasive_free"]],
      sae_total = e[["sae"]], sae_free = e[["sae_free"]], proc_deaths = e[["proc_deaths"]],
      dose_sv = e[["pet"]] * p$msv_pet / 1000, radiation_cancers = e[["rad_ca"]],
      missed_investigated = e[["missed_inv"]], qaly_total = l[["qaly"]] - s$q_harm, qaly_harms = s$q_harm,
      daly_total = l[["daly"]] + s$d_harm, daly_harms = s$d_harm)[ADDITIVE]
  }
  sapply(stats::setNames(arms, arms), one)
}

cohort_rows <- function(data, population = POPULATION) data$cohort[data$cohort$population == population, ]

run_population <- function(p, data, settings, arms, cohort = cohort_rows(data)) {
  settings <- kunz_settings(p, settings); cohort <- kunz_cohort(cohort, settings)
  wk <- workup_summary(p)
  mats <- lapply(seq_len(nrow(cohort)), function(i) row_outcomes(cohort[i, ], data, p, settings, arms, wk = wk))
  list(cohort = cohort, rows = mats, total = Reduce(`+`, Map(`*`, mats, cohort$share)))
}

# Galleri-specific derived metrics (ceakit adds inc_cost, inc_effect, icer, nmb)
derive_galleri <- function(A, settings) {
  per <- function(x, y) ifelse(y > 1e-15, x / y, NA)
  r <- A[, 1]; not_ref <- seq_len(ncol(A)) != 1
  d_av <- r["daly_total"] - A["daly_total", ]; ca_av <- r["ca_deaths"] - A["ca_deaths", ]
  rbind(mean_days = per(A["days_x_inv", ], A["investigated", ]),
        mean_days_fp = per(A["days_free", ], A["free_investigated", ]),
        overdiagnosis_share = per(A["overdiagnosed", ], A["screen_detected", ]),
        cost_tx = A["cost_total", ] - A["cost_tests", ],
        late_averted = r["late_lifetime"] - A["late_lifetime", ],
        ca_deaths_averted = ca_av, ly_gained = A["ly_total", ] - r["ly_total"],
        qaly_net = A["qaly_total", ] - r["qaly_total"], daly_averted = d_av,
        cost_per_daly = ifelse(not_ref, (A["cost_total", ] - r["cost_total"]) / d_av, NA),
        nn_per_ca_death_averted = ifelse(not_ref, per(1, ca_av), NA),
        investigated_per_found = per(A["investigated", ], A["screen_detected", ]),
        fp_per_found = per(A["free_investigated", ], A["screen_detected", ]),
        invasive_per_found = per(A["invasive_total", ], A["screen_detected", ]))
}

# ---- Subgroup grid: start age x sex x race/ethnicity (deterministic) --------
subgroup_cohort <- function(age, sex, race, settings) data.frame(
  population = POPULATION, race_eth = race, sex = sex, start_age = age,
  risk_group = c("kunzmann_high", "kunzmann_low"), share = c(kunz_share(settings, sex), 1 - kunz_share(settings, sex)),
  kunzmann_scoped = c(1, 0), stringsAsFactors = FALSE)

subgroup_grid <- function(p, data, settings, arms, cohort_n = 1e5) {
  settings <- kunz_settings(p, settings)
  g <- expand.grid(start_age = settings$subgroups$ages, sex = settings$subgroups$sexes,
                   race_eth = settings$subgroups$races, stringsAsFactors = FALSE)
  purrr::map_dfr(seq_len(nrow(g)), function(i) {
    C <- subgroup_cohort(g$start_age[i], g$sex[i], g$race_eth[i], settings)
    A <- run_population(p, data, settings, arms, C)$total
    m <- rbind(A, derive_galleri(A, settings)); ref <- m[, 1]
    tibble::tibble(start_age = g$start_age[i], sex = g$sex[i], race_eth = g$race_eth[i],
                   race_label = unname(RACE_LABELS[g$race_eth[i]]), arm = colnames(m),
                   cancers_found = m["screen_detected", ] * cohort_n,
                   false_positives = m["free_investigated", ] * cohort_n,
                   tests = m["galleri_tests", ] * cohort_n, true_pos = m["screen_tp", ] * cohort_n,
                   false_pos = m["screen_fp", ] * cohort_n, true_neg = m["screen_tn", ] * cohort_n,
                   false_neg = m["screen_fn", ] * cohort_n, overdiagnosed = m["overdiagnosed", ] * cohort_n,
                   fp_qalys_lost = m["qaly_fp", ] * cohort_n, fp_dalys = m["daly_fp", ] * cohort_n,
                   cancer_deaths = m["ca_deaths", ] * cohort_n,
                   cancer_deaths_averted = m["ca_deaths_averted", ] * cohort_n,
                   net_qalys = m["qaly_net", ] * cohort_n, net_dalys_averted = m["daly_averted", ] * cohort_n,
                   cost_per_person = m["cost_total", ], qaly_per_person = m["qaly_total", ],
                   inc_cost = m["cost_total", ] - ref[["cost_total"]],
                   inc_qaly = m["qaly_total", ] - ref[["qaly_total"]],
                   icer_vs_usual_care = ifelse(colnames(m) == arms[1], NA,
                                               (m["cost_total", ] - ref[["cost_total"]]) /
                                                 (m["qaly_total", ] - ref[["qaly_total"]])),
                   nmb_vs_usual_care = settings$wtp * (m["qaly_total", ] - ref[["qaly_total"]]) -
                     (m["cost_total", ] - ref[["cost_total"]]))
  })
}

# ---- Galleri accuracy against SEER-driven risk, by attained age, sex, race ---
# Sensitivity and false-positive rate are test properties (same for everyone; lower
# after round 1). How many true and false positives, and so PPV and NPV, follow from
# each group's SEER-calibrated preclinical upper GI cancer and other-cancer risk and
# its life table. One cohort per sex x race/ethnicity, screened yearly from start_age.
accuracy_table <- function(p, data, settings, start_age = 50, cohort_n = 1e5) {
  wk <- workup_summary(p)
  g <- expand.grid(sex = settings$subgroups$sexes, race_eth = settings$subgroups$races, stringsAsFactors = FALSE)
  purrr::map_dfr(seq_len(nrow(g)), function(i) {
    C <- subgroup_cohort(start_age, g$sex[i], g$race_eth[i], settings); tot <- NULL
    for (j in seq_len(nrow(C))) {
      row <- C[j, ]; prep <- prepare_row(row, data, p, settings)
      s <- simulate_arm("Galleri", row, prep$curve, prep$cal, prep$p, settings, wk, prep$V, prep$q)
      m <- s$ev_t[, c("tests", "tp", "fp", "fp_other_ca", "free_inv", "tn", "fn", "found_ugi", "overdx_ugi")] * row$share
      tot <- if (is.null(tot)) m else tot + m
    }
    age <- start_age + seq_len(nrow(tot)) - 1; tot <- tot[age <= settings$screen_stop_age, , drop = FALSE]
    age <- age[age <= settings$screen_stop_age]; d <- as.data.frame(tot) * cohort_n
    div <- function(a, b) ifelse(b > 0, a / b, NA_real_)
    tibble::tibble(sex = g$sex[i], race_eth = g$race_eth[i], race_label = unname(RACE_LABELS[g$race_eth[i]]), age = age,
      tests = d$tests, tp = d$tp, fp = d$fp, fp_other_cancer = d$fp_other_ca, fp_cancer_free = d$free_inv, tn = d$tn, fn = d$fn,
      prevalence_per_100k_tested = div(d$tp + d$fn, d$tests) * 1e5,
      sensitivity = div(d$tp, d$tp + d$fn), specificity = div(d$tn, d$tn + d$fp),
      ppv_upper_gi = div(d$tp, d$tp + d$fp), ppv_any_cancer = div(d$tp + d$fp_other_ca, d$tp + d$fp), npv = div(d$tn, d$tn + d$fn),
      cancer_free_per_cancer_found = div(d$free_inv, d$tp + d$fp_other_ca))
  })
}

# ---- Calibration check: model vs SEER incidence by age ----------------------
# For each subgroup row of the cohort: SEER target, modelled symptomatic
# incidence under usual care (should match), and total incidence (symptomatic
# + screen-detected) in each screening arm.
calibration_table <- function(p, data, settings, arms = ARMS, cohort = cohort_rows(data)) {
  wk <- workup_summary(p)
  purrr::map_dfr(seq_len(nrow(cohort)), function(i) {
    row <- cohort[i, ]; prep <- prepare_row(row, data, p, settings); curve <- prep$curve
    target <- target_ugi(curve, row, settings)
    purrr::map_dfr(arms, function(a) {
      s <- simulate_arm(a, row, curve, prep$cal, prep$p, settings, wk, prep$V)
      tr <- s$trace[seq_len(nrow(curve)), , drop = FALSE]
      alive <- rowSums(tr[, !STATES %in% DEAD, drop = FALSE])
      pp <- 1 - exp(-1 / prep$p$sojourn_years)
      flows <- t(vapply(seq_len(nrow(curve)), function(t) {
        P <- transition_matrix(curve$mort_other[t], prep$cal$onset[t], prep$p, prep$cal$pl[t])
        tr[t, "pre"] * c(P["pre", "early"] + P["pre", "late"], P["pre", "late"])
      }, numeric(2)))
      tibble::tibble(race_eth = row$race_eth, sex = row$sex, risk_group = row$risk_group, share = row$share,
                     arm = a, age = curve$age_from, seer = target, onset = prep$cal$onset,
                     model_clinical = flows[, 1] / alive, clipped = prep$cal$clipped,
                     seer_late_share = prep$cal$pl, model_late_share = flows[, 2] / pmax(flows[, 1], 1e-300))
    })
  })
}

# Cohort-weighted trace by arm
trace_table <- function(p, data, settings, arms) {
  C <- cohort_rows(data); wk <- workup_summary(p)
  purrr::map_dfr(seq_len(nrow(C)), function(i) {
    row <- C[i, ]; prep <- prepare_row(row, data, p, settings)
    purrr::map_dfr(arms, function(a) {
      tr <- simulate_arm(a, row, prep$curve, prep$cal, prep$p, settings, wk, prep$V)$trace
      tibble::as_tibble(tr) |> dplyr::mutate(cycle = seq_len(nrow(tr)) - 1, arm = a, share = row$share)
    })
  }) |>
    tidyr::pivot_longer(dplyr::all_of(STATES), names_to = "state", values_to = "prop") |>
    dplyr::mutate(state = sub("^pre_m$", "pre", sub("^lead_", "", state))) |>
    dplyr::group_by(arm, cycle, state) |> dplyr::summarise(prop = sum(prop * share), .groups = "drop") |>
    dplyr::mutate(state = factor(state, setdiff(STATES, c("pre_m", "lead_early", "lead_late")),
                                 c("Well", "Undiagnosed cancer", "Early-stage dx", "Late-stage dx",
                                   "Cancer death", "Other-cause death")))
}
