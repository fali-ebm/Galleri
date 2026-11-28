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

# pre / pre_ol: preclinical upper GI cancer other than EAC (squamous esophageal and gastric;
# all upper GI when miscan_eac = 0), early (stage I-II) / late (III-IV), the same two-stage
# process as EAC (miscan_rates with sojourn_years and the SEER late share of this pool).
# pre_m: preclinical cancer already missed by Galleri that stays undetectable by
# Galleri (scenario p_fn_persist > 0); otherwise it behaves exactly like pre.
# pre_ee / pre_el: preclinical esophageal adenocarcinoma (EAC), early / late stage, and
# lead_ee / lead_el its screen-detected lead-time states. They are used only in the MISCAN-ESO
# scenario (miscan_eac = 1, see miscan_rates()); otherwise no one enters them.
# early_e / late_e: diagnosed EAC with esophagus-specific survival (m_early_eac, m_late_eac).
# pre_f / pre_fe / pre_fl: cancers that start during the year and cause symptoms before the
# next yearly test (other upper GI; EAC early / late), so no test ever sees them
# (within_year_onset = 1). In continuous time a share of new cancers does this; without these
# states every cancer would meet at least one test.
STATES <- c("well", "pre", "pre_ol", "pre_m", "pre_ee", "pre_el", "pre_f", "pre_fe", "pre_fl", "lead_early", "lead_late",
            "lead_ee", "lead_el", "early", "late", "early_e", "late_e", "dead_cancer", "dead_other")
DX_EARLY <- c("early", "early_e"); DX_LATE <- c("late", "late_e")
PRE <- c("pre", "pre_ol", "pre_m", "pre_ee", "pre_el", "pre_f", "pre_fe", "pre_fl")
DEAD <- c("dead_cancer", "dead_other")
CHANNELS <- c("cost", "qaly", "daly", "ly", "clin_dx", "clin_late", "clin_dx_eac", "clin_late_eac", "ca_death", "oc_death", "symp_perf")
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

# MISCAN-ESO scenario: esophageal adenocarcinoma, and the rest of upper GI cancer.
target_eac <- function(curve, row, settings)
  curve$inc_eso * curve$eac_share * kunz_rr(row$risk_group, settings, row$sex)
target_other <- function(curve, row, settings) curve$inc_eso * (1 - curve$eac_share) + curve$inc_stom
late_other <- function(curve, p) {
  eso <- curve$inc_eso * (1 - curve$eac_share); stom <- curve$inc_stom
  pl <- (eso * curve$late_eso + stom * curve$late_stom) / (eso + stom)
  or <- p$or_late_clin %||% 1
  pl * or / (1 - pl + pl * or)
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

# ---- MISCAN-ESO natural history for esophageal adenocarcinoma (scenario) -----------
# In MISCAN-ESO (CISNET; Kong 2014) preclinical EAC moves through stages before symptoms.
# Here: preclinical early (PE) -> symptoms, early stage (rate hE), or -> preclinical late (PL)
# (rate lam); PL -> symptoms, late stage (rate hL). Exponential durations, as in MISCAN. The
# three rates follow from three targets:
#   S  mean sojourn, onset to symptoms (MISCAN-ESO: men 4.9 y, women 4.5 y)
#        S = 1 / (hE + lam) + L / hL
#   L  share late (regional/distant) at symptomatic diagnosis (SEER esophagus, subgroup)
#        L = lam / (hE + lam)
#   r  mean time in PE / mean time in PL (eac_dwell_ratio; 1 = equal, assumption)
#        r = hL / (hE + lam)
# Solution: hL = (r + L) / S, a = hE + lam = hL / r, hE = a (1 - L), lam = a L.
# A screen-detected EAC is early if found in PE and late if found in PL, so the stage of
# screen-detected EAC comes out of the natural history instead of p_late_screen.
miscan_rates <- function(S, L, r = 1) {
  stopifnot(S > 0, L > 0, L < 1, r > 0)
  hL <- (r + L) / S; a <- hL / r
  list(a = a, hE = a * (1 - L), lam = a * L, hL = hL, S = S, L = L, r = r)
}
# One-year transition probabilities of the PE/PL process (no other-cause death).
miscan_year <- function(m) {
  a <- m$a; b <- m$hL; ea <- exp(-a); eb <- exp(-b)
  pe_pl <- if (abs(a - b) < 1e-10) m$lam * ea else m$lam * (eb - ea) / (a - b)   # in PL at year end
  c(pe_pe = ea, pe_pl = pe_pl, pe_early = m$hE / a * (1 - ea),
    pe_late = m$lam / a * (1 - ea) - pe_pl, pl_pl = eb, pl_late = 1 - eb)
}
# Where a cancer that starts at a uniformly distributed time within the year is at the end
# of that year (within_year_onset = 1), averaging over the remaining time u ~ U(0, 1).
# g(x) = E[exp(-x u)] = (1 - exp(-x)) / x.
#   other upper GI (presentation rate h = 1 / sojourn): presents before year end with d = 1 - g(h)
#   EAC: still PE g(a); in PL lam (g(hL) - g(a)) / (a - hL); presented early hE / a (1 - g(a));
#        presented late lam / a (1 - g(a)) - (in PL)
onset_fractions <- function(p, eac, oth = eac) {
  one <- function(m) {
    if (!isTRUE((p$within_year_onset %||% 0) == 1)) return(c(pe = 1, pl = 0, fe = 0, fl = 0))
    g <- function(x) if (x < 1e-10) 1 - x / 2 else (1 - exp(-x)) / x
    a <- m$a; b <- m$hL
    pl <- if (abs(a - b) < 1e-8) m$lam * (1 - exp(-a) * (1 + a)) / a^2 else m$lam * (g(b) - g(a)) / (a - b)
    c(pe = g(a), pl = pl, fe = m$hE / a * (1 - g(a)), fl = m$lam / a * (1 - g(a)) - pl)
  }
  e <- one(eac); o <- one(oth); names(o) <- paste0("o", names(o)); c(e, o)
}

# EAC targets for one subgroup row: sex-specific sojourn S, dwell ratio r and the
# incidence-weighted SEER late share L. Other-cause death removes more people on the longer
# path to late-stage symptoms, so the late share seen at diagnosis is below lam / (hE + lam);
# prepare_row() therefore solves the death-free share so that the calibrated usual-care
# cohort reproduces L exactly (eac_late_share()).
eac_targets <- function(p, curve, row, settings) {
  S <- if (identical(row$sex, "female")) p$sojourn_eac_female %||% 4.5 else p$sojourn_eac_male %||% 4.9
  w <- curve$inc_eso * curve$eac_share; L <- sum(w * curve$late_eso) / sum(w)
  or <- p$or_late_clin %||% 1
  list(S = S, r = p$eac_dwell_ratio %||% 1, L = L * or / (1 - L + L * or))
}
eac_rates <- function(p, curve, row, settings) {           # rates with the death-free share solved
  tg <- eac_targets(p, curve, row, settings)
  calib_eac(curve, row, settings, p, tg)$eac
}
# Calibrate both cancers and solve the death-free late share Ls so that the lifetime late
# share of symptomatic EAC in the usual-care cohort equals the SEER target.
# Both two-stage processes (the non-EAC pool, and EAC when miscan_eac = 1) are calibrated
# together: onset by age to SEER incidence (calibrate) and the death-free late share Ls of each
# so that the lifetime late share at symptomatic diagnosis in the usual-care cohort equals the
# SEER target (incidence-weighted over ages). Fixed-point iteration on Ls (the map is close
# to the identity, so a few steps reach 1e-9).
calib_all <- function(curve, row, settings, p) {
  mis <- miscan_on(p); r <- p$eac_dwell_ratio %||% 1
  if (mis) { T_o <- target_other(curve, row, settings); w <- T_o; Lc <- late_other(curve, p)
             tg <- eac_targets(p, curve, row, settings); T_e <- target_eac(curve, row, settings)
  } else { T_o <- target_ugi(curve, row, settings); w <- T_o; Lc <- target_late(curve, row, settings, p); tg <- NULL; T_e <- NULL }
  L_o <- sum(w * Lc) / sum(w)
  Lso <- L_o; Lse <- if (mis) tg$L else NA
  clampL <- function(v) min(0.9999, max(0.0005, v))
  for (it in 1:40) {
    cal <- calibrate(curve, T_o, p, miscan_rates(p$sojourn_years, Lso, r), T_e, if (mis) miscan_rates(tg$S, Lse, tg$r))
    go <- L_o - cal$o_late / cal$o_clin; ge <- if (mis) tg$L - cal$eac_late / cal$eac_clin else 0
    if (abs(go) < 1e-9 && abs(ge) < 1e-9) break
    Lso <- clampL(Lso + go); if (mis) Lse <- clampL(Lse + ge)
  }
  cal$oth$L_seen <- L_o; if (mis) cal$eac$L_seen <- tg$L
  cal$fit_iter <- it
  cal
}
calib_eac <- function(curve, row, settings, p, tg = NULL) calib_all(curve, row, settings, p)
miscan_on <- function(p) isTRUE((p$miscan_eac %||% 0) == 1)

# ---- Galleri sensitivity by stage and cancer type -----------------------------------
# Galleri reports 12-month EPISODE sensitivity: true positives / (true positives + cancers
# diagnosed within 12 months of a negative test). The model needs the probability that a
# test is positive for a cancer that is preclinical at the test, by preclinical stage
# (early = stage I-II, late = III-IV) and type. Each cancer type is a two-stage
# preclinical process (miscan_rates: sojourn S, SEER late share L, dwell ratio r):
#   share of prevalent cases in early  piE = (1/a) / (1/a + (lam/a)/hL), late piL = 1 - piE
#   presents within 12 months if missed: from early qE = 1 - [exp(-a) + lam (exp(-hL) - exp(-a)) / (a - hL)]
#                                        from late  qL = 1 - exp(-hL)
#   episode sensitivity = (sE piE + sL piL) / (sE piE + sL piL + (1 - sE) piE qE + (1 - sL) piL qL)
# Unknowns: late-stage sensitivity for esophagus and stomach and one early/late ratio rho
# (sE = rho sL). Targets (GRAIL FDA panel summary 2026, PMA Galleri):
#   sens_ep_eso  esophagus episode sensitivity (pooled PATHFINDER 2 + NHS-Galleri)
#   sens_ep_stom stomach episode sensitivity
#   screen_early_share  share of Galleri-detected upper GI cancers at stage I-II (Figure 12)
# The episode sensitivities are multiplied by cso_acc (the share given the right CSO, here
# stomach/esophagus, which is what sends a person to EGD). Reference population: US, all
# races, both sexes, ages 50-79 (the trials' ages), prevalence = SEER incidence x sojourn.
stage_props <- function(m) {
  a <- m$a; b <- m$hL
  piE <- (1 / a) / (1 / a + (m$lam / a) / b)
  still <- exp(-a) + (if (abs(a - b) < 1e-10) m$lam * exp(-a) else m$lam * (exp(-b) - exp(-a)) / (a - b))
  c(piE = piE, piL = 1 - piE, qE = 1 - still, qL = 1 - exp(-b))
}
sens_ref_types <- function(p, data) {
  B <- data$baseline_risks; B <- B[B$race_eth == "all" & B$age_from >= 50 & B$age_from <= 79, ]
  inc <- c(eac = sum(B$inc_eso * B$eac_share), escc = sum(B$inc_eso * (1 - B$eac_share)), gas = sum(B$inc_stom))
  L <- c(eac = sum(B$inc_eso * B$eac_share * B$late_eso) / inc[["eac"]],
         escc = sum(B$inc_eso * (1 - B$eac_share) * B$late_eso) / inc[["escc"]],
         gas = sum(B$inc_stom * B$late_stom) / inc[["gas"]])
  wm <- sum(B$inc_eso[B$sex == "male"] * B$eac_share[B$sex == "male"]) / inc[["eac"]]
  S <- c(eac = wm * (p$sojourn_eac_male %||% 4.9) + (1 - wm) * (p$sojourn_eac_female %||% 4.5),
         escc = p$sojourn_years, gas = p$sojourn_years)
  r <- p$eac_dwell_ratio %||% 1
  lapply(setNames(names(inc), names(inc)), function(k) list(prev = inc[[k]] * S[[k]], st = stage_props(miscan_rates(S[[k]], L[[k]], r))))
}
galleri_stage_sens <- function(p, data) {
  ty <- sens_ref_types(p, data); acc <- p$cso_acc %||% 1
  tg <- c(eso = p$sens_ep_eso * acc, stom = p$sens_ep_stom * acc); early_t <- p$screen_early_share
  groups <- list(eso = c("eac", "escc"), stom = "gas")
  eps <- function(g, sL, rho) {
    num <- den <- 0
    for (k in groups[[g]]) { s <- ty[[k]]$st; w <- ty[[k]]$prev; sE <- min(1, rho * sL)
      d <- sE * s[["piE"]] + sL * s[["piL"]]
      num <- num + w * d; den <- den + w * (d + (1 - sE) * s[["piE"]] * s[["qE"]] + (1 - sL) * s[["piL"]] * s[["qL"]]) }
    num / den
  }
  solve_sL <- function(g, rho) {
    if (tg[[g]] <= 0) return(0)                                    # a test that detects nothing
    hi <- min(1, 1 / max(rho, 1e-9))
    f <- function(sL) eps(g, sL, rho) - tg[[g]]
    if (f(hi) <= 0) return(hi)                                   # cannot reach the target: cap
    stats::uniroot(f, c(1e-9, hi), tol = 1e-12)$root
  }
  early <- function(rho) {
    num <- den <- 0
    for (g in names(groups)) { sL <- solve_sL(g, rho)
      for (k in groups[[g]]) { s <- ty[[k]]$st; w <- ty[[k]]$prev; sE <- min(1, rho * sL)
        num <- num + w * sE * s[["piE"]]; den <- den + w * (sE * s[["piE"]] + sL * s[["piL"]]) } }
    if (den > 0) num / den else 0
  }
  f <- function(rho) early(rho) - early_t
  rho <- if (early_t < 0 || all(tg <= 0)) 1 else if (f(1e-6) >= 0) 1e-6 else if (f(20) <= 0) 20 else stats::uniroot(f, c(1e-6, 20), tol = 1e-12)$root   # early_t < 0: equal sensitivity
  sL <- c(eso = solve_sL("eso", rho), stom = solve_sL("stom", rho))
  list(rho = rho, sE = setNames(pmin(1, rho * sL), names(sL)), sL = sL, early = early(rho),
       eps = c(eso = eps("eso", sL[["eso"]], rho), stom = eps("stom", sL[["stom"]], rho)), target = tg)
}
# Per-cycle sensitivity and stage-at-detection for one subgroup. The cohort keeps the
# non-EAC upper GI cancers (squamous esophageal and gastric; also EAC when miscan_eac = 0) in
# one preclinical pool, so the pool's detection probability and late share at detection
# are mixed over types (weights: SEER incidence x sojourn of this subgroup and age), each type at
# its steady-state stage mix. EAC with miscan_eac = 1 uses its own early / late states.
#   s_pool   P(Galleri E/S positive | pool cancer)       pl_det  late share of Galleri-found pool cancers
#   pl_egd   late share of pool cancers found by an EGD (no Galleri; e.g. Kunzmann arm)
#   s_ee, s_el  EAC early / late; w_*  type shares of the pool (prevalence; Galleri-found early / late)
stage_sens_table <- function(curve, row, settings, p, gs) {
  mis <- miscan_on(p); r <- p$eac_dwell_ratio %||% 1; rr <- kunz_rr(row$risk_group, settings, row$sex)
  S_eac <- if (identical(row$sex, "female")) p$sojourn_eac_female %||% 4.5 else p$sojourn_eac_male %||% 4.9
  n <- nrow(curve); out <- vector("list", n)
  for (t in seq_len(n)) {
    ty <- list(escc = list(inc = curve$inc_eso[t] * (1 - curve$eac_share[t]), L = curve$late_eso[t], S = p$sojourn_years, g = "eso"),
               gas = list(inc = curve$inc_stom[t], L = curve$late_stom[t], S = p$sojourn_years, g = "stom"))
    if (!mis) ty$eac <- list(inc = curve$inc_eso[t] * curve$eac_share[t] * rr, L = curve$late_eso[t], S = p$sojourn_years, g = "eso")
    det_e <- det_l <- pre_ <- prl <- numeric(0)
    for (k in names(ty)) { z <- ty[[k]]; st <- stage_props(miscan_rates(z$S, min(max(z$L, 0.001), 0.999), r))
      w <- z$inc * z$S; pre_[k] <- w * st[["piE"]]; prl[k] <- w * st[["piL"]]
      det_e[k] <- pre_[k] * gs$sE[[z$g]]; det_l[k] <- prl[k] * gs$sL[[z$g]] }
    out[[t]] <- c(s_pool_e = sum(det_e) / sum(pre_), s_pool_l = sum(det_l) / sum(prl),
                  s_ee = gs$sE[["eso"]], s_el = gs$sL[["eso"]],
                  setNames(pre_ / sum(pre_), paste0("w_pre_", names(pre_))), setNames(prl / sum(prl), paste0("w_prl_", names(prl))),
                  setNames(det_e / max(sum(det_e), 1e-300), paste0("w_dete_", names(det_e))),
                  setNames(det_l / max(sum(det_l), 1e-300), paste0("w_detl_", names(det_l))))
  }
  as.data.frame(do.call(rbind, lapply(out, function(v) v[order(names(v))])))
}

transition_matrix <- function(mo, onset, p, oth, onset_e = 0, eac = NULL) {
  s <- 1 - mo
  n <- length(STATES)
  P <- matrix(0, n, n, dimnames = list(STATES, STATES))
  P[, "dead_other"] <- mo
  if (is.null(eac)) eac <- oth                                   # unused unless onset_e > 0
  fr <- onset_fractions(p, eac, oth)
  P["well", c("well", "pre", "pre_ol", "pre_f", "pre_ee", "pre_el", "pre_fe", "pre_fl")] <-
    s * c(1 - onset - onset_e, onset * fr[["ope"]], onset * fr[["opl"]], onset * (fr[["ofe"]] + fr[["ofl"]]),
          onset_e * fr[["pe"]], onset_e * fr[["pl"]], onset_e * fr[["fe"]], onset_e * fr[["fl"]])
  fq <- fr[["ofe"]] + fr[["ofl"]]; ffl <- if (fq > 0) fr[["ofl"]] / fq else oth$L
  P["pre_f", c("early", "late")] <- s * c(1 - ffl, ffl)
  P["pre_fe", "early_e"] <- s; P["pre_fl", "late_e"] <- s
  attr(P, "frac") <- fr
  yo <- miscan_year(oth)
  P["pre", c("pre", "pre_ol", "early", "late")] <- s * yo[c("pe_pe", "pe_pl", "pe_early", "pe_late")]
  P["pre_ol", c("pre_ol", "late")] <- s * yo[c("pl_pl", "pl_late")]
  # pre_m (scenario p_fn_persist): a cancer Galleri keeps missing; it still progresses, so it presents
  # at the pool's own stage mix (death-free late share L) over its mean sojourn S
  hm <- 1 - exp(-1 / oth$S)
  P["pre_m", c("pre_m", "early", "late")] <- s * c(1 - hm, hm * (1 - oth$L), hm * oth$L)
  qo <- exp(-1 / oth$S); qol <- exp(-oth$hL)                     # lead time: mean remaining sojourn
  P["lead_early", c("lead_early", "early")] <- s * c(qo, 1 - qo)
  P["lead_late", c("lead_late", "late")] <- s * c(qol, 1 - qol)
  y <- miscan_year(eac)
  P["pre_ee", c("pre_ee", "pre_el", "early_e", "late_e")] <- s * y[c("pe_pe", "pe_pl", "pe_early", "pe_late")]
  P["pre_el", c("pre_el", "late_e")] <- s * y[c("pl_pl", "pl_late")]
  # lead time of a screen-detected EAC: exponential with the mean remaining sojourn
  # (S from PE, 1 / hL from PL), so the would-be symptom time keeps MISCAN's mean
  ql <- exp(-1 / eac$S); qll <- exp(-eac$hL)
  P["lead_ee", c("lead_ee", "early_e")] <- s * c(ql, 1 - ql)
  P["lead_el", c("lead_el", "late_e")] <- s * c(qll, 1 - qll)
  me <- p$m_early_eac %||% p$m_early; ml <- p$m_late_eac %||% p$m_late
  P["early_e", c("early_e", "dead_cancer")] <- s * c(1 - me, me)
  P["late_e", c("late_e", "dead_cancer")] <- s * c(1 - ml, ml)
  P["early", c("early", "dead_cancer")] <- s * c(1 - p$m_early, p$m_early)
  P["late", c("late", "dead_cancer")] <- s * c(1 - p$m_late, p$m_late)
  P[DEAD, ] <- 0; P["dead_cancer", "dead_cancer"] <- 1; P["dead_other", "dead_other"] <- 1
  P
}

# Solve yearly onset so that modelled symptomatic incidence per living person
# equals the SEER target at every age (closed form, one year at a time).
# Returns the onset vector, the starting state (preclinical prevalence implied
# by the target) and the fit.
#
# MISCAN-ESO scenario (target_e and eac given): esophageal adenocarcinoma has its own
# target and its own two-stage preclinical process; 'target' is then the rest of upper GI
# cancer (gastric and other esophageal). Both onsets are solved the same way: the EAC onset
# puts enough people in PE that next year's EAC presentations (from PE and PL) per living
# person equal target_e. Without target_e the result is identical to the original model.
calibrate <- function(curve, target, p, oth, target_e = NULL, eac = NULL) {
  n <- nrow(curve); s <- 1 - curve$mort_other
  use_e <- !is.null(target_e)
  x <- setNames(numeric(length(STATES)), STATES)
  yo <- miscan_year(oth); oE <- yo[["pe_early"]] + yo[["pe_late"]]; oL <- yo[["pl_late"]]; ko <- oth$lam / oth$hL
  x["pre"] <- min(0.5, target[1] / (s[1] * (oE + ko * oL))); x["pre_ol"] <- ko * x[["pre"]]
  if (use_e) {
    y <- miscan_year(eac); cE <- y[["pe_early"]] + y[["pe_late"]]; cL <- y[["pl_late"]]
    k <- eac$lam / eac$hL                       # steady state: PL / PE = lam / hL
    x["pre_ee"] <- min(0.2, target_e[1] / (s[1] * (cE + k * cL))); x["pre_el"] <- k * x[["pre_ee"]]
  }
  x["well"] <- 1 - x[["pre"]] - x[["pre_ol"]] - x[["pre_ee"]] - x[["pre_el"]]
  x0 <- x; onset <- numeric(n); onset_e <- numeric(n); clipped <- 0; eac_clin <- 0; eac_late <- 0; o_clin <- 0; o_late <- 0
  prev <- matrix(0, n, 4, dimnames = list(NULL, c("o", "ol", "ee", "el")))
  for (t in seq_len(n)) {
    un <- sum(x[c("well", PRE)]); prev[t, ] <- c(x[["pre"]], x[["pre_ol"]], x[["pre_ee"]], x[["pre_el"]]) / un
    if (t < n) {
      P0 <- transition_matrix(curve$mort_other[t], 0, p, oth, 0, eac)
      nxt <- drop(x %*% P0)                       # next state with no new onset
      alive_next <- sum(nxt[!STATES %in% DEAD])
      # next year's presentations: pre (rate pp) plus this year's onsets that present within the year
      fr <- onset_fractions(p, eac %||% oth, oth); W <- max(1e-12, nxt[["well"]])
      o <- (target[t + 1] * alive_next / s[t + 1] - oE * nxt[["pre"]] - oL * nxt[["pre_ol"]]) /
           (W * (oE * fr[["ope"]] + oL * fr[["opl"]] + fr[["ofe"]] + fr[["ofl"]]))
      if (o < 0 || o > 1) clipped <- clipped + 1
      onset[t] <- min(1, max(0, o))
      if (use_e) {
        oe <- (target_e[t + 1] * alive_next / s[t + 1] - cE * nxt[["pre_ee"]] - cL * nxt[["pre_el"]]) /
              (W * (cE * fr[["pe"]] + cL * fr[["pl"]] + fr[["fe"]] + fr[["fl"]]))
        if (oe < 0 || oe > 1) clipped <- clipped + 1
        onset_e[t] <- min(1 - onset[t], max(0, oe))
      }
    } else { onset[t] <- onset[max(1, t - 1)]; onset_e[t] <- onset_e[max(1, t - 1)] }
    Pt <- transition_matrix(curve$mort_other[t], onset[t], p, oth, onset_e[t], eac)
    op <- c("pre", "pre_ol", "pre_f")
    o_late <- o_late + sum(x[op] * Pt[op, "late"]); o_clin <- o_clin + sum(x[op] * (Pt[op, "early"] + Pt[op, "late"]))
    if (use_e) {
      ep <- c("pre_ee", "pre_el", "pre_fe", "pre_fl")
      eac_late <- eac_late + sum(x[ep] * Pt[ep, "late_e"])
      eac_clin <- eac_clin + sum(x[ep] * (Pt[ep, "early_e"] + Pt[ep, "late_e"]))
    }
    x <- drop(x %*% Pt)
  }
  list(onset = onset, onset_e = onset_e, x0 = x0, clipped = clipped, oth = oth, eac = eac,
       eac_clin = eac_clin, eac_late = eac_late, o_clin = o_clin, o_late = o_late,
       prev_o = prev[, "o"], prev_ol = prev[, "ol"], prev_ee = prev[, "ee"], prev_el = prev[, "el"])   # usual-care prevalence, start of cycle
}
# transition matrix of cycle t for a calibrated row
tm_at <- function(curve, cal, p, t)
  transition_matrix(curve$mort_other[t], cal$onset[t], p, cal$oth, cal$onset_e[t] %||% 0, cal$eac)

# ---- Rewards per cycle (as before; onset now comes from calibration) --------
cycle_rewards <- function(P, curve, t, p) {
  alive <- setNames(ceakit::alive_fraction(P, DEAD), STATES)
  live <- !STATES %in% DEAD
  R <- matrix(0, length(STATES), length(CHANNELS), dimnames = list(STATES, CHANNELS))
  u <- c(well = 1, pre = 1, pre_ol = 1, pre_m = 1, pre_ee = 1, pre_el = 1, pre_f = 1, pre_fe = 1, pre_fl = 1, lead_early = p$u_mult_early, lead_late = p$u_mult_late,
         lead_ee = p$u_mult_early, lead_el = p$u_mult_late,
         early = p$u_mult_early, late = p$u_mult_late, early_e = p$u_mult_early, late_e = p$u_mult_late,
         dead_cancer = 0, dead_other = 0)
  dw <- c(well = 0, pre = 0, pre_ol = 0, pre_m = 0, pre_ee = 0, pre_el = 0, pre_f = 0, pre_fe = 0, pre_fl = 0, lead_early = p$dw_early, lead_late = p$dw_late,
          lead_ee = p$dw_early, lead_el = p$dw_late,
          early = p$dw_early, late = p$dw_late, early_e = p$dw_early, late_e = p$dw_late, dead_cancer = 0, dead_other = 0)
  R[, "ly"] <- alive
  R[, "qaly"] <- alive * curve$utility_well[t] * u[STATES]
  R["lead_early", "cost"] <- alive[["lead_early"]] * p$c_fu_early
  R["lead_late", "cost"]  <- alive[["lead_late"]] * p$c_fu_late
  R["lead_ee", "cost"] <- alive[["lead_ee"]] * p$c_fu_early
  R["lead_el", "cost"] <- alive[["lead_el"]] * p$c_fu_late
  R["early", "cost"] <- alive[["early"]] * p$c_fu_early + P["early", "dead_cancer"] * p$c_eol
  R["late", "cost"]  <- alive[["late"]] * p$c_fu_late  + P["late", "dead_cancer"]  * p$c_eol
  R["early_e", "cost"] <- alive[["early_e"]] * p$c_fu_early + P["early_e", "dead_cancer"] * p$c_eol
  R["late_e", "cost"]  <- alive[["late_e"]] * p$c_fu_late  + P["late_e", "dead_cancer"]  * p$c_eol
  R[, "daly"] <- alive * dw[STATES] + P[, "dead_cancer"] * curve$daly_ref_le[t] * live
  # symptomatic diagnosis: workup (incl. an EGD that can perforate) and first-year treatment
  perf <- p$p_perf_egd %||% 0; harm_perf <- perf * (1 - (p$p_death_perf %||% 0)) * (p$d_sae %||% 0)
  for (st in PRE) {
    pe <- sum(P[st, DX_EARLY]); pl_ <- sum(P[st, DX_LATE]); dx <- pe + pl_
    R[st, "cost"] <- pe * (p$c_symp_dx + p$c_tx_first_early) +
                     pl_ * (p$c_symp_dx + p$c_tx_first_late) + dx * perf * (p$c_perf %||% 0)
    R[st, "qaly"] <- R[st, "qaly"] - dx * harm_perf
    R[st, "daly"] <- R[st, "daly"] + dx * harm_perf
    R[st, "clin_dx"] <- dx
    R[st, "clin_dx_eac"] <- sum(P[st, c("early_e", "late_e")]); R[st, "clin_late_eac"] <- P[st, "late_e"]
    R[st, "clin_late"] <- pl_
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

# Averages over the rest of the workup: for people with another cancer (paths weighted
# by where cancers are, probability x (found + missed)) and for people with no cancer outside the upper GI tract
# (paths weighted by their cancer-free share), plus where the latter's
# pancreatobiliary history ends up.
REST_COLS <- c("n_ent", "n_pet", "n_eus", "n_samp", "n_tb", "n_lap", "n_ww", "n_incid", "cost", "days", "sae", "rad", "invasive", "dfp")
rest_summary <- function(p, pb) {
  d <- rest_tree(p, pb)
  wc <- d$prob * (d$found + d$missed); a <- colSums(d[REST_COLS] * wc) / sum(wc)   # people with another cancer: paths where cancers are
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
# Returns a matrix (cycles x pre, ee, el): the same probability for a cancer in pre, in
# preclinical early EAC and in preclinical late EAC (MISCAN-ESO scenario).
overdx_from <- function(Ps) {
  n <- length(Ps); q <- matrix(1, n + 1, 4, dimnames = list(NULL, c("pre", "ol", "ee", "el")))
  for (t in n:1) {
    P <- Ps[[t]]
    q[t, "ol"] <- P["pre_ol", "dead_other"] + P["pre_ol", "pre_ol"] * q[t + 1, "ol"]
    q[t, "pre"] <- P["pre", "dead_other"] + P["pre", "pre"] * q[t + 1, "pre"] + P["pre", "pre_ol"] * q[t + 1, "ol"]
    q[t, "el"] <- P["pre_el", "dead_other"] + P["pre_el", "pre_el"] * q[t + 1, "el"]
    q[t, "ee"] <- P["pre_ee", "dead_other"] + P["pre_ee", "pre_ee"] * q[t + 1, "ee"] + P["pre_ee", "pre_el"] * q[t + 1, "el"]
  }
  q[seq_len(n), , drop = FALSE]
}
overdx_prob <- function(curve, cal, p) overdx_from(lapply(seq_len(nrow(curve)), function(t) tm_at(curve, cal, p, t)))

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
# ee / el: preclinical esophageal adenocarcinoma, early / late (MISCAN-ESO scenario only).
HCOLS <- c("well", "pre", "ol", "iv", "pm", "ee", "el", "f", "fe", "fl")   # ol: pool, preclinical late   # f*: present before the next test
# k is tracked exactly up to KMAX false positives (KMAX = "KMAX or more"); with persistent false
# positives (fp_kappa small) a few people are positive almost every year, so the cap must be high.
KMAX <- 12
tag_grid <- function(H) expand.grid(k = 0:KMAX, h = 0:H, b = 1:3)
tag_id <- function(T, k, h, b) which(T$k == k & T$h == h & T$b == b)

# One screening round of the Galleri arm, all history tags at once (vectorised).
tag_idx <- function(k, h, b, H) 1 + k + (KMAX + 1) * h + (KMAX + 1) * (H + 1) * (b - 1)   # row of tag_grid(H)
acc_rows <- function(z, j, v) { r <- rowsum(v, j); i <- as.integer(rownames(r)); z[i] <- z[i] + r[, 1]; z }
# EGD rule (egd_rule): 0 = fixed hold, no EGD for egd_hold_years after a negative EGD;
# 1 = Bayesian hold: after a negative EGD, a repeat positive gets an EGD only if the posterior
# probability of an upper GI cancer, P(C | +, negative EGD y years ago, k earlier false
# positives), reaches the decision threshold (egd_threshold_at). The prior recovers from the
# EGD as new cancers arise: pi(y) = sum_c prev_c (1 - sens_egd exp(-y / S_c)), c = other upper
# GI (sojourn S_o) and EAC (S_e); P(+ | C) = Galleri sensitivity this round; P(+ | no C) =
# the tag's false-positive propensity plus positives from other-site cancers. Tags keep the
# years since the last negative EGD for egd_memory_years (h counts down; y = H - h + 1).
egd_posterior <- function(y, theta, t, nprev, curve, p, ctx) {
  rep_m <- if (nprev > 0) p$sens_g_rep_ratio %||% 1 else 1
  se <- p$sens_egd
  po_y <- ctx$prev_o * (1 - se * exp(-y / ctx$S_o)); pe_y <- ctx$prev_e * (1 - se * exp(-y / ctx$S_e))
  tp <- (po_y * ctx$s_u + pe_y * ctx$s_e) * rep_m; pi_y <- po_y + pe_y
  po <- curve$inc_other[t] * p$sens_g_other * rep_m
  out <- tp / (tp + (1 - pi_y) * (theta + po)); out[!is.finite(out)] <- 0; out
}
# Decision threshold for the EGD (QALYs, from the model's own values): do the EGD when
# P(C) x G >= (1 - P(C)) x L, i.e. P(C) >= L / (L + G), with
#   G = QALYs gained by finding a preclinical cancer now (lead-time state vs staying undiagnosed),
#       prevalence-weighted over other upper GI and EAC early / late;
#   L = QALYs an EGD costs a person without cancer: d_proc + perforation (non-fatal d_sae;
#       fatal: the person's remaining QALYs).
# egd_threshold >= 0 fixes the threshold instead.
egd_threshold_at <- function(p, cal, V, t) {
  if ((p$egd_threshold %||% -1) >= 0) return(p$egd_threshold)
  v <- V[[t]][, "qaly"]
  wo <- cal$prev_o[t] * cal$gs$s_pool_e[t]; wol <- cal$prev_ol[t] * cal$gs$s_pool_l[t]
  go <- (wo * (v[["lead_early"]] - v[["pre"]]) + wol * (v[["lead_late"]] - v[["pre_ol"]])) / max(1e-15, wo + wol)
  # EAC: weighted by the Galleri-detected stage mix (prevalence x stage sensitivity)
  we <- cal$prev_ee[t] * cal$gs$s_ee[t]; wl <- cal$prev_el[t] * cal$gs$s_el[t]
  ge <- (we * (v[["lead_ee"]] - v[["pre_ee"]]) + wl * (v[["lead_el"]] - v[["pre_el"]])) / max(1e-15, we + wl)
  po <- cal$prev_o[t] + cal$prev_ol[t]; pe <- cal$prev_ee[t] + cal$prev_el[t]
  G <- if (po + pe > 0) (po * go + pe * ge) / (po + pe) else go
  L <- p$d_proc + p$p_perf_egd * ((1 - p$p_death_perf) * p$d_sae + p$p_death_perf * v[["well"]])
  L / (L + max(G, 1e-9))
}
# Context for one Galleri round: usual-care prevalence (for the EGD rule), sojourns and the
# stage-specific sensitivities of this cycle.
egd_ctx <- function(cal, p, t, thr = NA) {
  g <- cal$gs[t, ]; pe <- cal$prev_ee[t]; pl <- cal$prev_el[t]; po <- cal$prev_o[t]; pol <- cal$prev_ol[t]
  list(prev_o = po + pol, prev_e = pe + pl, S_o = p$sojourn_years, S_e = (cal$eac %||% list(S = p$sojourn_years))$S,
       s_ue = g$s_pool_e, s_ul = g$s_pool_l, s_ee = g$s_ee, s_el = g$s_el,
       s_u = if (po + pol > 0) (po * g$s_pool_e + pol * g$s_pool_l) / (po + pol) else g$s_pool_e,
       s_e = if (pe + pl > 0) (pe * g$s_ee + pl * g$s_el) / (pe + pl) else g$s_ee, thr = thr)
}
# Map a vector named by cancer type (eac / escc / gas, any subset) onto all three.
type_vec <- function(v) { out <- c(eac = 0, escc = 0, gas = 0); nm <- sub("^w_(dete|detl|pre|prl)_", "", names(v)); out[nm] <- out[nm] + v; out }
# Split the pool's symptomatic early / late diagnoses of cycle t by type (SEER incidence x stage share).
clin_type_w <- function(curve, row, settings, p, t) {
  rr <- kunz_rr(row$risk_group, settings, row$sex)
  inc <- c(eac = if (miscan_on(p)) 0 else curve$inc_eso[t] * curve$eac_share[t] * rr,
           escc = curve$inc_eso[t] * (1 - curve$eac_share[t]), gas = curve$inc_stom[t])
  L <- c(eac = curve$late_eso[t], escc = curve$late_eso[t], gas = curve$late_stom[t])
  list(early = inc * (1 - L) / sum(inc * (1 - L)), late = inc * L / sum(inc * L))
}
hold_years <- function(p) if ((p$egd_rule %||% 0) == 1) max(1, round(p$egd_memory_years %||% 6)) else max(0, round(p$egd_hold_years %||% 0))

galleri_round <- function(X, T, t, nprev, curve, p, settings, wk, ctx = NULL) {
  H <- max(T$h); interval <- max(1, settings$screen_interval); n <- nrow(T)
  rep_m <- if (nprev > 0) p$sens_g_rep_ratio %||% 1 else 1
  m <- p$fpr_g_es; kap <- p$fp_kappa %||% 1e6
  sE <- 1 - p$p_perf_egd * p$p_death_perf; fi <- p$f_invisible %||% 0; se_egd <- p$sens_egd
  reach_pet <- 1 - p$p_ent_route; ppu <- p$p_pet_ugi %||% 0
  Fm <- do.call(rbind, lapply(wk$rest, `[[`, "free"))[T$b, , drop = FALSE]   # per tag: cancer-free averages
  Am <- do.call(rbind, lapply(wk$rest, `[[`, "all"))[T$b, , drop = FALSE]    # per tag: other-cancer averages
  det <- vapply(wk$rest, `[[`, 0, "det_o")[T$b]
  TO <- do.call(rbind, lapply(wk$rest, `[[`, "to"))[T$b, , drop = FALSE]
  W <- X[, "well"]; U <- X[, "pre"]; Ul <- X[, "ol"]; I <- X[, "iv"]; Mm <- X[, "pm"]; Ee <- X[, "ee"]; El <- X[, "el"]
  k <- T$k; h <- T$h; b <- T$b
  theta <- (m * kap + k) / (kap + nprev)                           # posterior mean false-positive propensity
  pool_o <- W * curve$inc_other[t] * (if (t == 1) p$sojourn_years else interval)
  a_u <- U * ctx$s_ue * rep_m; a_i <- I * ctx$s_ue * rep_m              # pool, early (I-II): stage I-II sensitivity
  a_ul <- Ul * ctx$s_ul * rep_m                                         # pool, late (III-IV)
  a_o <- pool_o * p$sens_g_other * rep_m; a_f <- (W - pool_o) * theta
  a_ee <- Ee * ctx$s_ee * rep_m; a_el <- El * ctx$s_el * rep_m          # EAC: early / late stage sensitivity
  gp <- a_u + a_ul + a_i + a_o + a_f + a_ee + a_el
  egd <- if (is.null(ctx) || (p$egd_rule %||% 0) != 1) h == 0 else
    h == 0 | egd_posterior(H - h + 1, theta, t, nprev, curve, p, ctx) >= ctx$thr
  s_e <- ifelse(egd, sE, 1)
  fu1 <- ifelse(egd, a_u * se_egd, 0)
  mv <- ifelse(egd, a_u * (1 - se_egd) * (1 - fi), a_u)          # visible, missed (or EGD skipped)
  mi <- ifelse(egd, a_u * (1 - se_egd) * fi, 0) + a_i              # not yet visible
  fu2 <- mv * s_e * reach_pet * ppu                                # seen on PET -> targeted EGD
  # pool, late stage: misses stay visible (as EAC)
  ful1 <- ifelse(egd, a_ul * se_egd, 0); mvol <- a_ul - ful1; ful2 <- mvol * s_e * reach_pet * ppu
  # EAC: found at EGD by stage; every EAC miss stays visible (in its own stage) and can be seen on PET
  fe1 <- ifelse(egd, a_ee * se_egd, 0); fl1 <- ifelse(egd, a_el * se_egd, 0)
  mve <- a_ee - fe1; mvl <- a_el - fl1
  fe2 <- mve * s_e * reach_pet * ppu; fl2 <- mvl * s_e * reach_pet * ppu
  fo <- a_o * s_e * det
  n_egd <- ifelse(egd, gp, 0) + fu2 + ful2 + fe2 + fl2
  rest_free <- (a_f + mv + mi + mve + mvl + mvol) * s_e; rest_o <- a_o * s_e
  R <- function(col) sum(rest_free * Fm[, col] + rest_o * Am[, col])
  perf <- n_egd * p$p_perf_egd; sae_o_v <- rest_free * Fm[, "sae"] + rest_o * Am[, "sae"]
  pd <- perf * p$p_death_perf + sae_o_v * p$p_death_sae
  ff <- a_f * s_e; pf <- ifelse(egd, a_f * p$p_perf_egd, 0); af_e <- ifelse(egd, a_f, 0)
  FF <- function(col) sum(ff * Fm[, col])
  e <- setNames(as.list(numeric(length(EVENTS))), EVENTS)
  Fq <- X[, "f"] + X[, "fe"] + X[, "fl"]                            # start after this test: tested as cancer-free
  e$tests <- sum(W + U + Ul + I + Mm + Ee + El + Fq); e$investigated <- sum(gp); e$free_inv <- sum(a_f)
  fe <- fe1 + fe2; fl <- fl1 + fl2
  e$found_ugi <- sum(fu1 + fu2 + ful1 + ful2 + fe + fl); e$found_other <- sum(fo); e$found_pet_ugi <- sum(fu2 + ful2 + fe2 + fl2)
  e$missed_inv <- sum(a_u + a_ul + a_i + a_o + a_ee + a_el - fu1 - fu2 - ful1 - ful2 - fe - fl - fo); e$egd_skipped <- sum(gp[!egd])
  e$egd <- sum(n_egd); e$egd_pos <- sum(fu1 + ful1 + fe1 + fl1); e$egd_neg <- sum(gp[egd]) - sum(fu1 + ful1 + fe1 + fl1)
  for (kk in c("ent", "pet", "eus", "samp", "tb", "lap")) e[[kk]] <- R(paste0("n_", kk))
  e$surv <- R("n_ww"); e$incid <- R("n_incid")
  e$invasive <- sum(n_egd) + R("invasive"); e$perf <- sum(perf); e$perf_deaths <- sum(perf) * p$p_death_perf
  e$sae <- sum(perf) + sum(sae_o_v); e$rad_ca <- R("rad")
  e$days <- sum(gp[egd]) * p$t_egd + R("days")
  e$tp <- sum(a_u + a_ul + a_i + a_ee + a_el); e$fp <- sum(a_o + a_f); e$fp_other_ca <- sum(a_o); e$tn <- sum(W - a_o - a_f + Fq)   # Fq: cancer starts after this test
  e$fn <- sum(U + Ul + I - a_u - a_ul - a_i + Mm + Ee + El - a_ee - a_el)
  e$hidden_new <- sum(mi); e$missed_visible <- sum(mv - fu2 + mvol - ful2 + mve - fe2 + mvl - fl2); e$repeat_fp <- sum(a_f[k > 0])
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
  cost_t <- sum(W + U + Ul + I + Mm + Ee + El + Fq) * p$c_galleri * (p$galleri_cost_share %||% 1) + sum(n_egd) * p$c_egd + R("cost")
  # where everyone goes next round
  h_after <- pmax(h - 1, 0); h_pos <- ifelse(egd, H, h_after)
  stay <- tag_idx(k, h_after, b, H)
  Xn <- X * 0
  Xn[, "well"] <- acc_rows(Xn[, "well"], stay, W - a_o - a_f - pd)
  Xn[, "pre"]  <- acc_rows(Xn[, "pre"], stay, U - a_u - persist)
  Xn[, "iv"]   <- acc_rows(Xn[, "iv"], stay, I - a_i)
  Xn[, "ol"]   <- acc_rows(Xn[, "ol"], stay, Ul - a_ul)
  Xn[, "pm"]   <- acc_rows(Xn[, "pm"], stay, Mm + persist)
  for (cc in c("f", "fe", "fl")) Xn[, cc] <- acc_rows(Xn[, cc], stay, X[, cc])
  Xn[, "ee"]   <- acc_rows(Xn[, "ee"], stay, Ee - a_ee)
  Xn[, "el"]   <- acc_rows(Xn[, "el"], stay, El - a_el)
  for (b2 in 1:3) {                                                # positives not diagnosed: hold and history updated
    w2 <- TO[, b2]
    j_fp <- tag_idx(pmin(k + 1, KMAX), h_pos, b2, H); j <- tag_idx(k, h_pos, b2, H)
    Xn[, "well"] <- acc_rows(Xn[, "well"], j_fp, a_f * w2)
    Xn[, "well"] <- acc_rows(Xn[, "well"], j, a_o * w2)
    Xn[, "pre"]  <- acc_rows(Xn[, "pre"], j, (mv - fu2) * w2)
    Xn[, "iv"]   <- acc_rows(Xn[, "iv"], j, mi * w2)
    Xn[, "ol"]   <- acc_rows(Xn[, "ol"], j, (mvol - ful2) * w2)
    Xn[, "ee"]   <- acc_rows(Xn[, "ee"], j, (mve - fe2) * w2)
    Xn[, "el"]   <- acc_rows(Xn[, "el"], j, (mvl - fl2) * w2)
  }
  list(X = Xn, e = lapply(e, function(v) unname(v)), q_path = unname(q_path), fp_cost = unname(fp_cost),
       cost_t = unname(cost_t), sae_other = unname(sum(sae_o_v)), pd = unname(sum(pd)),
       fe = unname(sum(fe)), fl = unname(sum(fl)), fue = unname(sum(fu1 + fu2)), ful = unname(sum(ful1 + ful2)))
}

# Natural history for the tagged people (the same transitions as transition_matrix);
# returns the next tagged matrix and the flows into diagnosed and dead states.
tagged_step <- function(X, T, P, onset, p, onset_e = 0) {
  s <- 1 - P["well", "dead_other"]
  v <- 1 - exp(-1 / (p$t_visible %||% 1))
  # hr_ca_after_fp: a false positive marks people already at higher risk; it does not
  # cause cancer. Onset is raised in people with an earlier false positive and lowered
  # in the rest so that total onset stays at the calibrated (usual care) level.
  hr <- p$hr_ca_after_fp %||% 1; w1 <- sum(X[T$k > 0, "well"]); w0 <- sum(X[T$k == 0, "well"])
  r0 <- if (hr != 1 && w0 > 0) max(0, (w0 + w1 - hr * w1) / w0) else 1
  mult <- ifelse(T$k > 0, hr, r0)
  on <- pmin(1, onset * mult); on_e <- pmin(1 - on, onset_e * mult)
  fr <- attr(P, "frac")
  # pool: early (pre, iv) and late (ol); iv stays invisible while early and becomes visible at rate v
  E1 <- X[, "pre"] + X[, "iv"]
  c_early <- E1 * P["pre", "early"] + X[, "pm"] * P["pre_m", "early"] + X[, "f"] * P["pre_f", "early"]
  c_late <- E1 * P["pre", "late"] + X[, "pm"] * P["pre_m", "late"] + X[, "ol"] * P["pre_ol", "late"] + X[, "f"] * P["pre_f", "late"]
  Xn <- X * 0
  Xn[, "well"] <- X[, "well"] * s * (1 - on - on_e)
  Xn[, "pre"] <- X[, "well"] * on * P["well", "pre"] / max(onset, 1e-300) + X[, "pre"] * P["pre", "pre"] + X[, "iv"] * P["pre", "pre"] * v
  Xn[, "ol"] <- X[, "well"] * on * P["well", "pre_ol"] / max(onset, 1e-300) + E1 * P["pre", "pre_ol"] + X[, "ol"] * P["pre_ol", "pre_ol"]
  Xn[, "f"] <- X[, "well"] * on * P["well", "pre_f"] / max(onset, 1e-300)
  Xn[, "fe"] <- X[, "well"] * s * on_e * fr[["fe"]]; Xn[, "fl"] <- X[, "well"] * s * on_e * fr[["fl"]]
  Xn[, "iv"] <- X[, "iv"] * P["pre", "pre"] * (1 - v)
  Xn[, "pm"] <- X[, "pm"] * P["pre_m", "pre_m"]
  # esophageal adenocarcinoma, MISCAN-ESO scenario (P rows already include survival s)
  Xn[, "ee"] <- X[, "well"] * s * on_e * fr[["pe"]] + X[, "ee"] * P["pre_ee", "pre_ee"]
  Xn[, "el"] <- X[, "well"] * s * on_e * fr[["pl"]] + X[, "ee"] * P["pre_ee", "pre_el"] + X[, "el"] * P["pre_el", "pre_el"]
  e_early <- sum(X[, "ee"]) * P["pre_ee", "early_e"] + sum(X[, "fe"]) * s
  e_late <- sum(X[, "ee"]) * P["pre_ee", "late_e"] + sum(X[, "el"]) * P["pre_el", "late_e"] + sum(X[, "fl"]) * s
  list(X = Xn, early = unname(sum(c_early)), late = unname(sum(c_late)),
       early_e = unname(e_early), late_e = unname(e_late), dead_other = unname(sum(X) * (1 - s)))
}

agg_state <- function(X, y) { x <- y; x[["well"]] <- sum(X[, "well"]); x[["pre"]] <- sum(X[, "pre"]) + sum(X[, "iv"]); x[["pre_ol"]] <- sum(X[, "ol"])
  x[["pre_m"]] <- sum(X[, "pm"]); x[["pre_ee"]] <- sum(X[, "ee"]); x[["pre_el"]] <- sum(X[, "el"])
  x[["pre_f"]] <- sum(X[, "f"]); x[["pre_fe"]] <- sum(X[, "fe"]); x[["pre_fl"]] <- sum(X[, "fl"]); x }

simulate_arm <- function(arm, row, curve, cal, p, settings, wk, V, q_od = overdx_prob(curve, cal, p)) {
  n <- nrow(curve); pls_o <- p$p_late_screen; x <- cal$x0
  cost_fp <- 0; qaly_fp <- 0
  fpq <- c(pathway = 0, worry = 0, complications = 0, deaths = 0, radiation = 0)   # QALYs lost by cancer-free people worked up
  fpd <- fpq; names(fpd)[1] <- "anxiety"                                            # DALYs, same people
  ev <- setNames(numeric(length(EVENTS)), EVENTS)
  acc <- setNames(numeric(length(CHANNELS)), CHANNELS)
  cost0 <- 0; q_harm <- 0; d_harm <- 0
  st_scr <- c(pool_early = 0, pool_late = 0, eac_early = 0, eac_late = 0)       # screen-found upper GI by stage
  st_type <- matrix(0, 3, 4, dimnames = list(c("eac", "escc", "gas"), c("scr_early", "scr_late", "clin_early", "clin_late")))
  trace <- matrix(0, n + 1, length(STATES), dimnames = list(NULL, STATES)); trace[1, ] <- x
  ev_t <- matrix(0, n, length(EVENTS), dimnames = list(NULL, EVENTS))        # screening events by cycle
  acc_t <- matrix(0, n, length(CHANNELS), dimnames = list(NULL, CHANNELS))   # undiscounted channel flows by cycle
  interval <- max(1, settings$screen_interval)
  is_g <- arm == "Galleri"
  if (is_g) {                                                   # tagged people (history), all start untagged
    T <- tag_grid(hold_years(p))
    X <- matrix(0, nrow(T), length(HCOLS), dimnames = list(NULL, HCOLS))
    X[tag_idx(0, 0, 1, max(T$h)), c("well", "pre", "ol", "pm", "ee", "el", "f", "fe", "fl")] <-
      c(x[["well"]], x[["pre"]], x[["pre_ol"]], x[["pre_m"]], x[["pre_ee"]], x[["pre_el"]], x[["pre_f"]], x[["pre_fe"]], x[["pre_fl"]])
  }
  for (t in seq_len(n)) {
    age <- row$start_age + t - 1
    df <- disc_vec(settings, t)
    P <- tm_at(curve, cal, p, t)
    if (arm == "Kunzmann -> EGD" && t == 1) { ev["kunz"] <- 1; cost0 <- cost0 + df[["cost"]] * p$c_kunz }
    screen <- arm != "Usual care" && age <= settings$screen_stop_age && (t - 1) %% interval == 0 &&
              (is_g || (isTRUE(row$kunzmann_scoped == 1) &&
                        ev[["rounds"]] < (settings$kunz_egd_rounds %||% Inf)))   # Kunzmann: one-time EGD by default
    if (screen) {
      if (is_g) {
        W <- sum(X[, "well"])
        ctx <- egd_ctx(cal, p, t, thr = egd_threshold_at(p, cal, V, t))
        g <- galleri_round(X, T, t, ev[["rounds"]], curve, p, settings, wk, ctx)
        X <- g$X; e <- g$e; sae_other <- g$sae_other; fp_cost <- g$fp_cost; cost_t <- g$cost_t
        pd <- g$pd; fe <- g$fe; fl <- g$fl; fue <- g$fue; ful <- g$ful
      } else {
        W <- x[["well"]]; U <- x[["pre"]] + x[["pre_ol"]]; e <- setNames(as.list(numeric(length(EVENTS))), EVENTS)
        fue <- x[["pre"]] * p$sens_egd; ful <- x[["pre_ol"]] * p$sens_egd   # pool, early / late
        Ue <- x[["pre_ee"]]; Ul <- x[["pre_el"]]                    # EAC
        fe <- Ue * p$sens_egd; fl <- Ul * p$sens_egd
        egd <- W + U + Ue + Ul
        e$found_ugi <- fue + ful + fe + fl; e$missed_inv <- U + Ue + Ul - e$found_ugi
        e$investigated <- egd; e$free_inv <- W
        e$egd <- egd; e$invasive <- egd; e$invasive_free <- W
        e$perf <- egd * p$p_perf_egd; e$perf_deaths <- e$perf * p$p_death_perf; sae_other <- 0
        e$egd_free <- W; e$perf_free <- W * p$p_perf_egd; e$deaths_free <- e$perf_free * p$p_death_perf
        e$days_free <- W * p$t_egd
        fp_cost <- W * (p$c_egd + p$p_perf_egd * p$c_perf)
        e$sae <- e$perf; e$sae_free <- W * p$p_perf_egd; e$days <- egd * p$t_egd
        e$egd_pos <- e$found_ugi; e$egd_neg <- egd - e$found_ugi   # EGD is the screening test here
        e$tp <- e$found_ugi; e$fn <- U + Ue + Ul - e$found_ugi; e$tn <- W
        cost_t <- egd * p$c_egd
        pd <- min(W * 0.5, e$perf_deaths)
      }
      e$fn_symp <- e$fn * (P["pre", "early"] + P["pre", "late"] + P["pre", "pre_ol"] * P["pre_ol", "late"])   # false negatives presenting with symptoms within a year
      # found upper GI cancers keep the stage they were found in: pool early / late (fue / ful),
      # EAC early / late (fe / fl). Type split of the pool: stage_sens_table weights.
      fu <- fue + ful; fo <- e$found_other
      gsr <- cal$gs[t, ]
      wde <- unlist(gsr[grep(if (is_g) "^w_dete_" else "^w_pre_", names(gsr))]); wdl <- unlist(gsr[grep(if (is_g) "^w_detl_" else "^w_prl_", names(gsr))])
      st_scr <- c(st_scr[1:2] + c(fue, ful), st_scr[3:4] + c(fe, fl))
      st_type[, "scr_early"] <- st_type[, "scr_early"] + type_vec(fue * wde) + c(eac = fe, escc = 0, gas = 0)
      st_type[, "scr_late"] <- st_type[, "scr_late"] + type_vec(ful * wdl) + c(eac = fl, escc = 0, gas = 0)
      e$overdx_ugi <- fue * q_od[t, "pre"] + ful * q_od[t, "ol"] + fe * q_od[t, "ee"] + fl * q_od[t, "el"]
      e$overdx_other <- fo * q_od[t, "pre"]
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
      e$proc_deaths <- pd; e$late_screen <- ful + fo * pls_o + fl; e$rounds <- 1
      for (k in names(e)) { ev[k] <- ev[k] + e[[k]]; ev_t[t, k] <- e[[k]] }
      # screen-detected upper GI cancers into lead time; procedure deaths out of well
      if (!is_g) { x["pre"] <- x["pre"] - fue; x["pre_ol"] <- x["pre_ol"] - ful; x["pre_ee"] <- x["pre_ee"] - fe; x["pre_el"] <- x["pre_el"] - fl
                   x["well"] <- x["well"] - pd }
      x["lead_early"] <- x["lead_early"] + fue; x["lead_late"] <- x["lead_late"] + ful
      x["lead_ee"] <- x["lead_ee"] + fe; x["lead_el"] <- x["lead_el"] + fl
      x["dead_other"] <- x["dead_other"] + pd
      if (is_g) x <- agg_state(X, x)
      tx_first <- fue * p$c_tx_first_early + ful * p$c_tx_first_late +
                  fo * ((1 - pls_o) * p$c_tx_first_early + pls_o * p$c_tx_first_late) +
                  fe * p$c_tx_first_early + fl * p$c_tx_first_late
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
        v_det <- (1 - pls_o) * V[[t]]["lead_early", ] + pls_o * V[[t]]["lead_late", ]
        delta <- (v_det - V[[t]]["pre", ]) * df
        acc <- acc + fo * delta
      }
    }
    R <- cycle_rewards(P, curve, t, p)
    acc_t[t, ] <- drop(x %*% R)
    cw <- clin_type_w(curve, row, settings, p, t); a_ <- acc_t[t, ]
    pe_ <- (a_[["clin_dx"]] - a_[["clin_late"]]) - (a_[["clin_dx_eac"]] - a_[["clin_late_eac"]])
    pl_ <- a_[["clin_late"]] - a_[["clin_late_eac"]]
    st_type[, "clin_early"] <- st_type[, "clin_early"] + pe_ * cw$early + c(eac = a_[["clin_dx_eac"]] - a_[["clin_late_eac"]], escc = 0, gas = 0)
    st_type[, "clin_late"] <- st_type[, "clin_late"] + pl_ * cw$late + c(eac = a_[["clin_late_eac"]], escc = 0, gas = 0)
    acc <- acc + acc_t[t, ] * df
    if (is_g) {
      st <- tagged_step(X, T, P, cal$onset[t], p, cal$onset_e[t] %||% 0)
      y <- x; y[c("well", "pre", "pre_ol", "pre_m", "pre_ee", "pre_el", "pre_f", "pre_fe", "pre_fl")] <- 0
      y <- drop(y %*% P)
      y[["early"]] <- y[["early"]] + st$early; y[["late"]] <- y[["late"]] + st$late
      y[["early_e"]] <- y[["early_e"]] + st$early_e; y[["late_e"]] <- y[["late_e"]] + st$late_e
      y[["dead_other"]] <- y[["dead_other"]] + st$dead_other
      X <- st$X; x <- agg_state(X, y)
    } else x <- drop(x %*% P)
    trace[t + 1, ] <- x
  }
  list(ev = ev, acc = acc, cost0 = cost0, q_harm = q_harm, d_harm = d_harm, trace = trace,
       pre0 = cal$x0[["pre"]] + cal$x0[["pre_ee"]] + cal$x0[["pre_el"]], cost_fp = cost_fp, qaly_fp = qaly_fp, fpq = fpq, fpd = fpd, ev_t = ev_t, acc_t = acc_t,
       st_scr = st_scr, st_type = st_type)
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
  "missed_investigated", "qaly_total", "qaly_harms", "daly_total", "daly_harms",
  "stage12_ugi", "stage34_ugi", "stage12_eso", "stage34_eso", "stage12_eac", "stage34_eac",
  "stage12_escc", "stage34_escc", "stage12_gastric", "stage34_gastric",
  "screen_stage12_ugi", "screen_stage34_ugi", "screen_stage12_eso", "screen_stage34_eso")

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
  cal <- calib_all(curve, row, settings, p)
  cal$gs <- stage_sens_table(curve, row, settings, p, galleri_stage_sens(p, data))
  Ps <- lapply(seq_len(nrow(curve)), function(t) tm_at(curve, cal, p, t))
  Rs <- lapply(seq_len(nrow(curve)), function(t) cycle_rewards(Ps[[t]], curve, t, p))
  list(curve = curve, cal = cal, V = values_by_cycle(Ps, Rs, settings), q = overdx_from(Ps), p = p)
}

row_outcomes <- function(row, data, p, settings, arms, prep = prepare_row(row, data, p, settings),
                         wk = workup_summary(p)) {
  p <- prep$p
  one <- function(a) {
    s <- simulate_arm(a, row, prep$curve, prep$cal, p, settings, wk, prep$V, prep$q); e <- s$ev; l <- s$acc
    ST <- s$st_type; S12 <- ST[, "scr_early"] + ST[, "clin_early"]; S34 <- ST[, "scr_late"] + ST[, "clin_late"]
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
      daly_total = l[["daly"]] + s$d_harm, daly_harms = s$d_harm,
      # stage at diagnosis (screen-found + symptomatic), lifetime: stage I-II = early (endoscopically
      # resectable or stage II), III-IV = late. Types: eac, escc (squamous esophageal), gas (gastric)
      stage12_ugi = sum(S12), stage34_ugi = sum(S34), stage12_eso = sum(S12[c("eac", "escc")]),
      stage34_eso = sum(S34[c("eac", "escc")]), stage12_eac = S12[["eac"]], stage34_eac = S34[["eac"]],
      stage12_escc = S12[["escc"]], stage34_escc = S34[["escc"]], stage12_gastric = S12[["gas"]], stage34_gastric = S34[["gas"]],
      screen_stage12_ugi = sum(ST[, "scr_early"]), screen_stage34_ugi = sum(ST[, "scr_late"]),
      screen_stage12_eso = sum(ST[c("eac", "escc"), "scr_early"]), screen_stage34_eso = sum(ST[c("eac", "escc"), "scr_late"]))[ADDITIVE]
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
        # primary comparative endpoint: extra stage I-II diagnoses vs usual care (do nothing)
        stage12_gained_ugi = A["stage12_ugi", ] - r["stage12_ugi"], stage12_gained_eso = A["stage12_eso", ] - r["stage12_eso"],
        stage34_averted_ugi = r["stage34_ugi"] - A["stage34_ugi", ], stage34_averted_eso = r["stage34_eso"] - A["stage34_eso", ],
        stage12_share_ugi = per(A["stage12_ugi", ], A["stage12_ugi", ] + A["stage34_ugi", ]),
        stage12_share_eso = per(A["stage12_eso", ], A["stage12_eso", ] + A["stage34_eso", ]),
        screen_stage12_share = per(A["screen_stage12_ugi", ], A["screen_stage12_ugi", ] + A["screen_stage34_ugi", ]),
        cost_per_stage12_gained = ifelse(not_ref, per(A["cost_total", ] - r["cost_total"], A["stage12_ugi", ] - r["stage12_ugi"]), NA),
        cost_per_stage12_gained_eso = ifelse(not_ref, per(A["cost_total", ] - r["cost_total"], A["stage12_eso", ] - r["stage12_eso"]), NA),
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
        P <- tm_at(curve, prep$cal, prep$p, t)
        pr <- c("pre", "pre_ol", "pre_ee", "pre_el", "pre_f", "pre_fe", "pre_fl")
        c(sum(tr[t, pr] * rowSums(P[pr, c(DX_EARLY, DX_LATE)])), sum(tr[t, pr] * rowSums(P[pr, DX_LATE])))
      }, numeric(2)))
      tibble::tibble(race_eth = row$race_eth, sex = row$sex, risk_group = row$risk_group, share = row$share,
                     arm = a, age = curve$age_from, seer = target, onset = prep$cal$onset,
                     model_clinical = flows[, 1] / alive, clipped = prep$cal$clipped,
                     seer_late_share = target_late(curve, row, settings, prep$p), model_late_share = flows[, 2] / pmax(flows[, 1], 1e-300))
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
    dplyr::mutate(state = dplyr::recode(state, pre_ol = "pre", pre_m = "pre", pre_ee = "pre", pre_el = "pre", pre_f = "pre", pre_fe = "pre",
                                        pre_fl = "pre", lead_early = "early",
                                        lead_late = "late", lead_ee = "early", lead_el = "late", early_e = "early", late_e = "late")) |>
    dplyr::group_by(arm, cycle, state) |> dplyr::summarise(prop = sum(prop * share), .groups = "drop") |>
    dplyr::mutate(state = factor(state, c("well", "pre", "early", "late", "dead_cancer", "dead_other"),
                                 c("Well", "Undiagnosed cancer", "Early-stage dx", "Late-stage dx",
                                   "Cancer death", "Other-cause death")))
}
