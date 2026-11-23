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
#        Kunzmann -> EGD : people flagged by the score have an EGD;
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

# ---- Galleri +ve workup decision tree (flowchart) --------------------------
workup_edges <- function(p) data.frame(
  from = c("galleri", "egd", "egd", "egd_neg", "egd_neg", "ent", "ent", "pet", "pet", "pet",
           "pb_pos", "eus", "eus", "eus"),
  to   = c("egd", "egd_pos", "egd_neg", "ent", "pet", "ent_pos", "ent_neg", "pb_pos", "incid", "pet_neg",
           "eus", "eus_pos", "eus_neg", "no_eus"),
  prob = c(1, p$p_egd_pos, 1 - p$p_egd_pos, p$p_ent_route, 1 - p$p_ent_route, p$p_ent_pos, 1 - p$p_ent_pos,
           p$p_pet_pb, p$p_pet_incid, 1 - p$p_pet_pb - p$p_pet_incid,
           1, p$p_eus_pos, p$p_eus_neg, 1 - p$p_eus_pos - p$p_eus_neg))

WORKUP_NODES <- data.frame(   # procedures done at each node (summed along each path)
  id     = c("egd", "ent", "pet", "eus", "eus_pos", "eus_neg", "ent_neg", "no_eus", "incid", "pet_neg"),
  n_egd  = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  n_ent  = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
  n_pet  = c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
  n_eus  = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
  n_samp = c(0, 0, 0, 0, 1, 1, 0, 0, 0, 0),
  n_tb   = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1),
  n_incid = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0))   # follow-up of an incidental PET finding (cancer or benign)

workup_paths <- function(p) {
  d <- ceakit::tree_paths(workup_edges(p), WORKUP_NODES)
  found <- c(egd_pos = 1, ent_pos = 1, eus_pos = 1, incid = p$p_incid_cancer)
  site <- c(egd_pos = "Upper GI", ent_pos = "Head & neck", eus_pos = "Pancreatobiliary", incid = "Other (incidental)")
  d$found <- ifelse(d$leaf %in% names(found), found[d$leaf], 0)
  d$missed <- ifelse(d$leaf %in% c("ent_neg", "eus_neg", "no_eus", "pet_neg"), p$p_miss, 0)
  d$site <- unname(site[d$leaf])
  add_path_details(d, p)
}

add_path_details <- function(d, p) {
  d$free <- 1 - d$found - d$missed
  d$n_invasive <- d$n_egd + d$n_eus + d$n_samp
  # serious complications after the EGD (EGD perforation is handled separately in workup_summary)
  d$sae <- d$n_ent * p$sae_ent + d$n_eus * p$sae_eus + d$n_samp * p$sae_eus_samp
  d$rad_ca <- d$n_pet * p$msv_pet / 1000 * p$risk_per_sv
  d$days <- d$n_egd * p$t_egd + d$n_ent * p$t_ent + d$n_pet * p$t_pet + d$n_eus * p$t_eus + d$n_tb * p$t_tb
  d$cost_dx <- d$n_egd * p$c_egd + d$n_ent * p$c_ent + d$n_pet * p$c_pet + d$n_eus * p$c_eus +
               d$n_samp * p$c_eus_samp + d$n_tb * p$c_tb + d$n_incid * (p$c_incid_fu %||% 0)
  d
}

# Per-positive averages of the workup tree, plus the share of cancers among
# positives that the workup finds: found / (found + missed).
# Every Galleri +ve has the EGD. The EGD can perforate (p_perf_egd); a fatal
# perforation (p_death_perf) ends the workup, so every step after the EGD is
# scaled by the share surviving it, s_egd. Cancers among the very few fatal
# perforations are not removed from the cancer pool (procedure deaths are taken
# from the cancer-free well state), a negligible approximation.
workup_summary <- function(p) {
  d <- workup_paths(p); s <- function(x) sum(d$prob * x)
  perf <- p$p_perf_egd; perf_death <- perf * p$p_death_perf; s_egd <- 1 - perf_death
  yf <- s(d$found); ym <- s(d$missed)
  after <- function(x, egd_part) egd_part + s_egd * (s(x) - egd_part)   # EGD itself happens for everyone
  c(n_egd = 1, n_ent = s_egd * s(d$n_ent), n_pet = s_egd * s(d$n_pet), n_eus = s_egd * s(d$n_eus),
    n_samp = s_egd * s(d$n_samp), n_tb = s_egd * s(d$n_tb), n_incid = s_egd * s(d$n_incid),
    n_invasive = after(d$n_invasive, 1),
    perf = perf, perf_death = perf_death, sae_other = s_egd * s(d$sae), sae = perf + s_egd * s(d$sae),
    rad_ca = s_egd * s(d$rad_ca), days = after(d$days, p$t_egd), cost_dx = after(d$cost_dx, p$c_egd),
    egd_pos = p$p_egd_pos * s_egd, egd_neg = (1 - p$p_egd_pos) * s_egd,
    det = if (yf + ym > 0) yf / (yf + ym) else 0,
    free_summary(d, p, perf, perf_death, s_egd))
}

# Averages for a CANCER-FREE person worked up (false positive): the workup paths
# weighted by each path's cancer-free share, so false positives follow the paths
# they actually take (e.g. mostly to PET -ve / ENT -ve), not the average path.
free_summary <- function(d, p, perf, perf_death, s_egd) {
  w <- d$prob * d$free; wf <- sum(w)
  a <- function(x) if (wf > 0) sum(w * x) / wf else 0
  sae_o <- s_egd * a(d$sae)
  c(f_n_ent = s_egd * a(d$n_ent), f_n_pet = s_egd * a(d$n_pet), f_n_eus = s_egd * a(d$n_eus),
    f_n_samp = s_egd * a(d$n_samp), f_n_tb = s_egd * a(d$n_tb), f_n_incid = s_egd * a(d$n_incid),
    f_n_invasive = 1 + s_egd * (a(d$n_invasive) - 1), f_sae_other = sae_o, f_sae = perf + sae_o,
    f_deaths = perf_death + sae_o * p$p_death_sae, f_rad_ca = s_egd * a(d$rad_ca),
    f_days = p$t_egd + s_egd * (a(d$days) - p$t_egd), f_cost_dx = p$c_egd + s_egd * (a(d$cost_dx) - p$c_egd))
}

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
EVENTS <- c("kunz", "tests", "investigated", "egd", "ent", "pet", "eus", "samp", "tb", "invasive",
            "invasive_free", "free_inv", "sae", "sae_free", "proc_deaths", "rad_ca", "days", "found_ugi",
            "found_other", "missed_inv", "late_screen", "rounds",
            "tp", "fp", "fp_other_ca", "tn", "fn", "fn_symp", "egd_pos", "egd_neg", "perf", "perf_deaths",
            "overdx_ugi", "overdx_other", "egd_free", "ent_free", "pet_free", "eus_free", "samp_free",
            "perf_free", "deaths_free", "rad_free", "days_free", "incid", "incid_free", "fn_persist")

simulate_arm <- function(arm, row, curve, cal, p, settings, wk, V, q_od = overdx_prob(curve, cal, p)) {
  n <- nrow(curve); pls <- p$p_late_screen; x <- cal$x0
  cost_fp <- 0; qaly_fp <- 0
  ev <- setNames(numeric(length(EVENTS)), EVENTS)
  acc <- setNames(numeric(length(CHANNELS)), CHANNELS)
  cost0 <- 0; q_harm <- 0; d_harm <- 0
  trace <- matrix(0, n + 1, length(STATES), dimnames = list(NULL, STATES)); trace[1, ] <- x
  interval <- max(1, settings$screen_interval)
  for (t in seq_len(n)) {
    age <- row$start_age + t - 1
    df <- disc_vec(settings, t)
    P <- transition_matrix(curve$mort_other[t], cal$onset[t], p, cal$pl[t])
    if (arm == "Kunzmann -> EGD" && t == 1) { ev["kunz"] <- 1; cost0 <- cost0 + df[["cost"]] * p$c_kunz }
    screen <- arm != "Usual care" && age <= settings$screen_stop_age && (t - 1) %% interval == 0 &&
              (arm == "Galleri" || isTRUE(row$kunzmann_scoped == 1))
    if (screen) {
      W <- x[["well"]]; U <- x[["pre"]]; e <- list()
      if (arm == "Galleri") {
        Um <- x[["pre_m"]]                                     # missed before and still undetectable
        rep_m <- if (ev[["rounds"]] > 0) p$sens_g_rep_ratio %||% 1 else 1   # repeat (incidence) rounds
        pool_o <- W * curve$inc_other[t] * (if (t == 1) p$sojourn_years else interval)
        a_u <- U * p$sens_g_ugi * rep_m; a_o <- pool_o * p$sens_g_other * rep_m; a_f <- (W - pool_o) * p$fpr_g_es
        gp <- a_u + a_o + a_f
        e$found_ugi <- a_u * wk[["det"]]; e$found_other <- a_o * wk[["det"]]
        e$missed_inv <- (a_u + a_o) - e$found_ugi - e$found_other
        e$tests <- W + U + Um; e$investigated <- gp; e$free_inv <- a_f
        for (k in c("egd", "ent", "pet", "eus", "samp", "tb")) e[[k]] <- gp * wk[[paste0("n_", k)]]
        e$invasive <- gp * wk[["n_invasive"]]; e$invasive_free <- a_f * wk[["f_n_invasive"]]
        e$sae <- gp * wk[["sae"]]; e$sae_free <- a_f * wk[["f_sae"]]; e$rad_ca <- gp * wk[["rad_ca"]]
        # burden in cancer-free people worked up (false positives with no cancer)
        e$egd_free <- a_f; e$ent_free <- a_f * wk[["f_n_ent"]]; e$pet_free <- a_f * wk[["f_n_pet"]]
        e$eus_free <- a_f * wk[["f_n_eus"]]; e$samp_free <- a_f * wk[["f_n_samp"]]
        e$perf_free <- a_f * wk[["perf"]]; e$deaths_free <- a_f * wk[["f_deaths"]]
        e$rad_free <- a_f * wk[["f_rad_ca"]]; e$days_free <- a_f * wk[["f_days"]]
        fp_cost <- a_f * (wk[["f_cost_dx"]] + wk[["perf"]] * p$c_perf + wk[["f_sae_other"]] * p$c_sae)   # f_cost_dx includes incidental follow-up
        e$days <- gp * wk[["days"]]
        e$perf <- gp * wk[["perf"]]; e$perf_deaths <- gp * wk[["perf_death"]]; sae_other <- gp * wk[["sae_other"]]
        e$egd_pos <- gp * wk[["egd_pos"]]; e$egd_neg <- gp * wk[["egd_neg"]]
        # Galleri 2x2 against the target condition (preclinical upper GI cancer):
        # TP = cancer & +ve; FP = no upper GI cancer & +ve (another cancer or none);
        # FN = cancer & -ve (stays undiagnosed: symptoms later, or found at a later round); TN = rest
        e$tp <- a_u; e$fp <- a_o + a_f; e$fp_other_ca <- a_o; e$fn <- U - a_u + Um; e$tn <- W - a_o - a_f
        e$incid <- gp * wk[["n_incid"]]; e$incid_free <- a_f * wk[["f_n_incid"]]
        e$fn_persist <- (U - a_u) * (p$p_fn_persist %||% 0)    # newly missed cancers Galleri will keep missing
        cost_t <- (W + U + Um) * p$c_galleri * (p$galleri_cost_share %||% 1) + gp * wk[["cost_dx"]]
      } else {
        egd <- W + U
        e$found_ugi <- U * p$sens_egd; e$found_other <- 0; e$missed_inv <- U - e$found_ugi
        e$tests <- 0; e$investigated <- egd; e$free_inv <- W
        e$egd <- egd; e$invasive <- egd; e$invasive_free <- W
        e$perf <- egd * p$p_perf_egd; e$perf_deaths <- e$perf * p$p_death_perf; sae_other <- 0
        e$egd_free <- W; e$ent_free <- 0; e$pet_free <- 0; e$eus_free <- 0; e$samp_free <- 0
        e$perf_free <- W * p$p_perf_egd; e$deaths_free <- e$perf_free * p$p_death_perf
        e$rad_free <- 0; e$days_free <- W * p$t_egd
        fp_cost <- W * (p$c_egd + p$p_perf_egd * p$c_perf)
        e$sae <- e$perf; e$sae_free <- W * p$p_perf_egd; e$rad_ca <- 0; e$days <- egd * p$t_egd
        e$egd_pos <- e$found_ugi; e$egd_neg <- egd - e$found_ugi   # EGD is the screening test here
        e$tp <- e$found_ugi; e$fp <- 0; e$fp_other_ca <- 0; e$fn <- U - e$found_ugi; e$tn <- W
        e$incid <- 0; e$incid_free <- 0; e$fn_persist <- 0
        cost_t <- egd * p$c_egd
      }
      e$fn_symp <- e$fn * (P["pre", "early"] + P["pre", "late"])   # false negatives presenting with symptoms within a year
      fu <- e$found_ugi; fo <- e$found_other
      e$overdx_ugi <- fu * q_od[t]; e$overdx_other <- fo * q_od[t]
      pd <- min(W * 0.5, e$perf_deaths + sae_other * p$p_death_sae)
      # QALYs lost by cancer-free people worked up: procedures, non-fatal complications,
      # radiation cancers, and the remaining (healthy) life of those who die
      n_proc_free <- e$egd_free + e$ent_free + e$pet_free + e$eus_free + e$incid_free
      qh_free <- n_proc_free * p$d_proc + (e$sae_free - e$deaths_free) * p$d_sae + e$rad_free * p$loss_rad_ca
      qaly_fp <- qaly_fp + df[["qaly"]] * (qh_free + e$deaths_free * V[[t]]["well", "qaly"])
      cost_fp <- cost_fp + df[["cost"]] * fp_cost
      e$proc_deaths <- pd; e$late_screen <- (fu + fo) * pls; e$rounds <- 1
      for (k in names(e)) ev[k] <- ev[k] + e[[k]]
      # move screen-detected upper GI cancers into lead time; procedure deaths out of well
      x["pre"] <- x["pre"] - fu - e$fn_persist; x["pre_m"] <- x["pre_m"] + e$fn_persist
      x["lead_early"] <- x["lead_early"] + fu * (1 - pls); x["lead_late"] <- x["lead_late"] + fu * pls
      x["well"] <- x["well"] - pd; x["dead_other"] <- x["dead_other"] + pd
      tx_first <- (fu + fo) * ((1 - pls) * p$c_tx_first_early + pls * p$c_tx_first_late)
      cost0 <- cost0 + df[["cost"]] * (cost_t + tx_first + e$perf * p$c_perf + sae_other * p$c_sae)
      n_proc <- e$egd + (e$ent %||% 0) + (e$pet %||% 0) + (e$eus %||% 0) + e$incid
      qh <- n_proc * p$d_proc + (e$sae - pd) * p$d_sae + e$rad_ca * p$loss_rad_ca
      q_harm <- q_harm + df[["qaly"]] * qh
      d_harm <- d_harm + df[["daly"]] * (qh + pd * curve$daly_ref_le[t])
      # cancers elsewhere found now vs presenting with symptoms later
      if (fo > 0) {
        v_det <- (1 - pls) * V[[t]]["lead_early", ] + pls * V[[t]]["lead_late", ]
        delta <- (v_det - V[[t]]["pre", ]) * df
        acc <- acc + fo * delta
      }
    }
    R <- cycle_rewards(P, curve, t, p)
    acc <- acc + drop(x %*% R) * df
    x <- drop(x %*% P); trace[t + 1, ] <- x
  }
  list(ev = ev, acc = acc, cost0 = cost0, q_harm = q_harm, d_harm = d_harm, trace = trace,
       pre0 = cal$x0[["pre"]], cost_fp = cost_fp, qaly_fp = qaly_fp)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Outcomes per cohort row ------------------------------------------------
ADDITIVE <- c("total_cancers", "cancer_free", "kunz_scores", "galleri_tests", "screen_rounds",
  "screen_tp", "screen_fp", "screen_fp_other_cancer", "screen_tn", "screen_fn", "fn_symptomatic_1y",
  "egd_positive", "egd_negative_workup", "perforations", "perforation_deaths",
  "overdiagnosed", "overdiagnosed_ugi", "overdiagnosed_other", "egd_free", "ent_free", "pet_ct_free", "eus_free",
  "eus_sampling_free", "perforations_free", "proc_deaths_free", "dose_sv_free", "radiation_cancers_free",
  "days_free", "cost_fp", "qaly_fp", "incidental_followups", "incidental_followups_free",
  "perforations_symptomatic", "fn_persistent", "investigated",
  "egd", "ent", "pet_ct", "eus", "eus_sampling", "tumor_board", "invasive_total", "days_x_inv",
  "cost_tests", "cost_total", "screen_detected", "found_upper_gi", "found_other",
  "late_lifetime", "clin_dx_lifetime", "ca_deaths", "oc_deaths", "ly_total", "free_investigated",
  "invasive_free", "sae_total", "sae_free", "proc_deaths", "dose_sv", "radiation_cancers",
  "missed_investigated", "qaly_total", "qaly_harms", "daly_total", "daly_harms")

prepare_row <- function(row, data, p, settings) {
  curve <- baseline_for(data, row, settings)
  cal <- calibrate(curve, target_ugi(curve, row, settings), p, target_late(curve, row, settings, p))
  Ps <- lapply(seq_len(nrow(curve)), function(t) transition_matrix(curve$mort_other[t], cal$onset[t], p, cal$pl[t]))
  Rs <- lapply(seq_len(nrow(curve)), function(t) cycle_rewards(Ps[[t]], curve, t, p))
  q <- numeric(length(Ps) + 1); q[length(q)] <- 1                    # overdiagnosis probability (see overdx_prob)
  for (t in rev(seq_along(Ps))) q[t] <- Ps[[t]]["pre", "dead_other"] + Ps[[t]]["pre", "pre"] * q[t + 1]
  list(curve = curve, cal = cal, V = values_by_cycle(Ps, Rs, settings), q = q[seq_along(Ps)])
}

row_outcomes <- function(row, data, p, settings, arms, prep = prepare_row(row, data, p, settings)) {
  wk <- workup_summary(p)
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
      perforations_symptomatic = l[["symp_perf"]], fn_persistent = e[["fn_persist"]], investigated = e[["investigated"]], egd = e[["egd"]], ent = e[["ent"]],
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
  mats <- lapply(seq_len(nrow(cohort)), function(i) row_outcomes(cohort[i, ], data, p, settings, arms))
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
      s <- simulate_arm(a, row, curve, prep$cal, p, settings, wk, prep$V)
      tr <- s$trace[seq_len(nrow(curve)), , drop = FALSE]
      alive <- rowSums(tr[, !STATES %in% DEAD, drop = FALSE])
      pp <- 1 - exp(-1 / p$sojourn_years)
      flows <- t(vapply(seq_len(nrow(curve)), function(t) {
        P <- transition_matrix(curve$mort_other[t], prep$cal$onset[t], p, prep$cal$pl[t])
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
      tr <- simulate_arm(a, row, prep$curve, prep$cal, p, settings, wk, prep$V)$trace
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
