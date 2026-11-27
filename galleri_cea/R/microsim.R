# Individual-level microsimulation (Nov 2026)
#
# Natural history is simulated person by person with hesim's individual-level
# continuous-time state-transition model (IndivCtstm): well -> preclinical upper GI
# cancer -> clinical presentation (early or late), with other-cause death competing.
# Hazards are piecewise constant by year and come from the same calibrated onset,
# sojourn, stage split and life table as the cohort model.
#
# Screening is then overlaid on each person's own history: annual Galleri tests with a
# person-level false-positive propensity (beta(m kappa, (1 - m) kappa), drawn once), the
# EGD hold, cancers missed or not yet visible at EGD, the rest of the workup sampled path
# by path, pancreatobiliary history flags, perforation and other procedure deaths; the
# Kunzmann arm has its one-time EGD. All arms run on the same simulated people (common
# random numbers), so each person's incremental outcome is observed directly.
#
# Cancers at other sites found by the workup are not part of the natural history; as in
# the cohort model, each one found is credited with its expected value (found now vs
# presenting with symptoms later, from backward induction).

MS_STATES <- c("well", "pre", "clin_early", "clin_late", "dead_other")

# ---- Continuous-time calibration -------------------------------------------------
# The cohort model's yearly onset assumes death comes first within each year. In
# continuous time the same probabilities give more onset and presentation before death
# (most at old ages), so onset hazards are re-solved here, year by year, so that the
# expected presentations per living person equal the same SEER targets. Within a year
# hazards are constant: well -> pre (onset), pre -> clinical (1 / sojourn), other death.
ms_calibrate <- function(prep, settings, row) {
  curve <- prep$curve; p <- prep$p; cal <- prep$cal; nT <- nrow(curve)
  target <- target_ugi(curve, row, settings)
  lo <- -log(1 - curve$mort_other); lp <- 1 / p$sojourn_years
  lce <- -log(1 - p$m_early); lcl <- -log(1 - p$m_late)
  W <- cal$x0[["well"]]; P <- cal$x0[["pre"]]; De <- 0; Dl <- 0
  r_on <- numeric(nT); clipped <- 0
  g <- function(l) if (abs(l) < 1e-12) 1 else (1 - exp(-l)) / l       # mean of exp(-l s) over the year
  for (t in seq_len(nT)) {
    l2 <- lo[t] + lp; alive0 <- W + P + De + Dl
    pres <- function(lon) { l1 <- lo[t] + lon
      lp * (P * g(l2) + if (lon > 0) lon * W * (g(l1) - g(l2)) / (l2 - l1 + 1e-15) else 0) }
    goal <- target[t] * alive0
    if (pres(0) >= goal) { lon <- 0; clipped <- clipped + 1 }
    else lon <- uniroot(function(z) pres(z) - goal, c(0, 50), tol = 1e-14)$root
    r_on[t] <- lon; l1 <- lo[t] + lon
    Pnew <- P * exp(-l2) + lon * W * (exp(-l1) - exp(-l2)) / (l2 - l1 + 1e-15)
    pr <- pres(lon); pl <- cal$pl[t]
    De <- De * exp(-(lo[t] + lce)) + pr * (1 - pl) * exp(-(lo[t] + lce) / 2)
    Dl <- Dl * exp(-(lo[t] + lcl)) + pr * pl * exp(-(lo[t] + lcl) / 2)
    W <- W * exp(-l1); P <- Pnew
  }
  list(r_on = pmax(r_on, 1e-12), clipped = clipped)
}

# ---- Natural history with hesim -----------------------------------------------
ms_natural_history <- function(prep, n, seed, settings, row) {
  curve <- prep$curve; cal <- prep$cal; p <- prep$p; nT <- nrow(curve)
  r <- function(x) pmax(x, 1e-12)
  r_on <- ms_calibrate(prep, settings, row)$r_on; r_oc <- r(-log(1 - curve$mort_other))
  r_pr <- 1 / p$sojourn_years
  rates <- list(r_on, r_oc, r(r_pr * (1 - cal$pl)), r(r_pr * cal$pl), r_oc)
  tmat <- rbind(c(NA, 1, NA, NA, 2), c(NA, NA, 3, 4, 5), rep(NA, 5), rep(NA, 5), rep(NA, 5))
  dimnames(tmat) <- list(MS_STATES, MS_STATES)
  tm <- seq_len(nT) - 1
  mk <- function(v) hesim::params_surv(
    coefs = setNames(lapply(v, function(z) matrix(log(z), 1, 1, dimnames = list(NULL, "(Intercept)"))), paste0("rate", seq_along(v))),
    aux = list(time = tm), dist = "pwexp")
  pars <- do.call(hesim::params_surv_list, lapply(rates, mk))
  set.seed(seed)
  pre0 <- cal$x0[["pre"]] / (cal$x0[["well"]] + cal$x0[["pre"]])
  start <- ifelse(runif(n) < pre0, 2L, 1L)
  hd <- hesim::hesim_data(strategies = data.table::data.table(strategy_id = 1),
                          patients = data.table::data.table(patient_id = seq_len(n)))
  td <- hesim::expand(hd, by = c("strategies", "patients")); td[, `(Intercept)` := 1]
  tr <- hesim::create_IndivCtstmTrans(pars, input_data = td, trans_mat = tmat, clock = "forward", start_state = start)
  m <- hesim::IndivCtstm$new(trans_model = tr)
  m$sim_disease(max_t = nT)
  dp <- m$disprog_
  T_on <- ifelse(start == 2L, 0, Inf); T_clin <- rep(Inf, n); late <- rep(FALSE, n); T_oc <- rep(Inf, n)
  i <- dp$to == 2; T_on[dp$patient_id[i]] <- dp$time_stop[i]
  i <- dp$to %in% 3:4; T_clin[dp$patient_id[i]] <- dp$time_stop[i]; late[dp$patient_id[i]] <- dp$to[i] == 4
  i <- dp$to == 5; T_oc[dp$patient_id[i]] <- dp$time_stop[i]
  # after presentation: other-cause death (same yearly hazards, memoryless) and cancer death
  # (constant yearly hazard by stage); one uniform per person shared by both stages and arms
  H <- c(0, cumsum(r_oc))
  draw_after <- function(t0) {                         # inverse cumulative hazard from t0
    out <- rep(Inf, length(t0)); ok <- is.finite(t0) & t0 < nT
    h0 <- approx(0:nT, H, t0[ok], rule = 2)$y
    target <- h0 - log(runif(sum(ok)))
    out[ok] <- ifelse(target >= H[nT + 1], Inf, approx(H, 0:nT, target, ties = "ordered")$y)
    out
  }
  T_oc_post <- draw_after(T_clin)
  data.frame(id = seq_len(n), T_on = T_on, T_clin = T_clin, late = late, T_oc = T_oc, T_oc_post = T_oc_post,
             u_cd = runif(n))
}

# ---- Screening overlay ---------------------------------------------------------
ms_paths <- function(wk) lapply(wk$rest, function(r) {
  d <- r$paths
  d$w_free <- d$prob * d$free / sum(d$prob * d$free)
  d$w_ca <- d$prob * (d$found + d$missed) / sum(d$prob * (d$found + d$missed))
  d$p_found_ca <- ifelse(d$found + d$missed > 0, d$found / (d$found + d$missed), 0)
  d
})

ms_screen <- function(nh, arm, prep, row, settings, wk, seed) {
  p <- prep$p; curve <- prep$curve; V <- prep$V; q <- prep$q; nT <- nrow(curve); N <- nrow(nh)
  set.seed(seed)
  H <- max(0, round(p$egd_hold_years %||% 0)); kap <- p$fp_kappa %||% 1e6; m <- p$fpr_g_es
  theta <- if (kap >= 1e6) rep(m, N) else rbeta(N, m * kap, (1 - m) * kap)
  z <- function() numeric(N)
  O <- list(tests = z(), positives = z(), fp = z(), fp_other = z(), egd = z(), ent = z(), pet = z(), eus = z(),
            samp = z(), tb = z(), lap = z(), ww = z(), incid = z(), perf = z(), sae = z(), days = z(),
            hidden = z(), egd_skipped = z(), found_other = z(), worry = z(),
            cost = z(), qaly = z(), daly = z(), cost_fp = z(), qaly_fp = z(), daly_fp = z())
  dx_t <- rep(Inf, N); dx_late <- rep(NA, N); T_pd <- rep(Inf, N); pm <- rep(FALSE, N)
  iv_until <- rep(-Inf, N); hold <- integer(N); k_fp <- integer(N); pb <- rep(1L, N)
  PT <- ms_paths(wk); pls <- p$p_late_screen
  sE_d <- p$p_death_perf
  disc <- function(t, r) (1 + r)^-(t - 1)
  rounds <- 0
  last <- min(nT - 1, settings$screen_stop_age - row$start_age)
  taus <- if (arm == "Galleri") seq(0, last, by = max(1, settings$screen_interval)) else
          if (arm == "Kunzmann -> EGD" && isTRUE(row$kunzmann_scoped == 1)) 0 else numeric(0)
  if (arm == "Kunzmann -> EGD") O$cost <- O$cost + p$c_kunz
  for (tau in taus) {
    t <- tau + 1; dfc <- disc(t, settings$disc_cost); dfq <- disc(t, settings$disc_qaly); dfd <- disc(t, settings$disc_daly)
    alive <- nh$T_oc > tau & T_pd > tau & dx_t > tau & nh$T_clin > tau
    idx <- which(alive); if (!length(idx)) next
    pre <- nh$T_on[idx] <= tau
    if (arm == "Kunzmann -> EGD") {                   # one-time EGD for a score >= 8
      O$egd[idx] <- O$egd[idx] + 1; O$days[idx] <- O$days[idx] + p$t_egd
      found <- pre & runif(length(idx)) < p$sens_egd
      perf <- runif(length(idx)) < p$p_perf_egd; pd <- perf & runif(length(idx)) < sE_d
      O$perf[idx] <- O$perf[idx] + perf
      cst <- p$c_egd + perf * p$c_perf
      O$cost[idx] <- O$cost[idx] + dfc * cst
      ql <- p$d_proc + (perf & !pd) * p$d_sae
      O$qaly[idx] <- O$qaly[idx] - dfq * ql; O$daly[idx] <- O$daly[idx] + dfd * ql
      fr <- !pre
      O$cost_fp[idx[fr]] <- O$cost_fp[idx[fr]] + dfc * cst[fr]; O$qaly_fp[idx[fr]] <- O$qaly_fp[idx[fr]] + dfq * ql[fr]
      O$daly_fp[idx[fr]] <- O$daly_fp[idx[fr]] + dfd * ql[fr]
      T_pd[idx[pd]] <- tau
      fi <- idx[found & !pd]; dx_t[fi] <- tau; dx_late[fi] <- runif(length(fi)) < pls
      next
    }
    # ---- Galleri ----
    O$tests[idx] <- O$tests[idx] + 1
    O$cost[idx] <- O$cost[idx] + dfc * p$c_galleri * (p$galleri_cost_share %||% 1)
    rep_m <- if (rounds > 0) p$sens_g_rep_ratio %||% 1 else 1; rounds <- rounds + 1
    nI <- length(idx); u <- runif(nI)
    ugi <- pre & !pm[idx]
    pos_u <- ugi & u < p$sens_g_ugi * rep_m
    newpm <- ugi & !pos_u & runif(nI) < (p$p_fn_persist %||% 0); pm[idx[newpm]] <- TRUE
    pool <- curve$inc_other[t] * (if (t == 1) p$sojourn_years else 1)
    oth <- !pre & runif(nI) < pool
    pos_o <- oth & runif(nI) < p$sens_g_other * rep_m
    pos_f <- !pre & !oth & runif(nI) < theta[idx]
    pos <- pos_u | pos_o | pos_f
    h_now <- hold[idx]; hold[idx] <- pmax(h_now - 1L, 0L)                 # the hold counts down every round
    if (!any(pos)) next
    j <- idx[pos]; nj <- length(j); type <- ifelse(pos_u[pos], "u", ifelse(pos_o[pos], "o", "f"))
    egd <- h_now[pos] == 0
    O$positives[j] <- O$positives[j] + 1; O$fp[j] <- O$fp[j] + (type == "f"); O$fp_other[j] <- O$fp_other[j] + (type == "o")
    O$egd_skipped[j] <- O$egd_skipped[j] + !egd
    cst <- numeric(nj); qp <- numeric(nj); dp <- numeric(nj); dead <- rep(FALSE, nj); found <- rep(FALSE, nj)
    nproc <- numeric(nj); days <- numeric(nj); dfp_v <- numeric(nj)
    # step 1: EGD
    vis <- tau >= iv_until[j]
    if (any(egd)) {
      e <- which(egd)
      O$egd[j[e]] <- O$egd[j[e]] + 1; cst[e] <- cst[e] + p$c_egd; nproc[e] <- nproc[e] + 1; days[e] <- days[e] + p$t_egd
      perf <- runif(length(e)) < p$p_perf_egd; pdth <- perf & runif(length(e)) < sE_d
      O$perf[j[e]] <- O$perf[j[e]] + perf; O$sae[j[e]] <- O$sae[j[e]] + perf
      cst[e] <- cst[e] + perf * p$c_perf; qp[e] <- qp[e] + (perf & !pdth) * p$d_sae; dead[e] <- pdth
      ue <- e[type[e] == "u" & !pdth]
      if (length(ue)) {
        r <- runif(length(ue))
        f1 <- vis[ue] & r < p$sens_egd
        hid <- vis[ue] & !f1 & runif(length(ue)) < (p$f_invisible %||% 0)
        found[ue[f1]] <- TRUE
        iv_until[j[ue[hid]]] <- tau + rexp(sum(hid), 1 / (p$t_visible %||% 1))
        O$hidden[j[ue[hid]]] <- O$hidden[j[ue[hid]]] + 1
        O$hidden[j[ue[!vis[ue]]]] <- O$hidden[j[ue[!vis[ue]]]] + 1           # still not visible
      }
      neg <- e[!found[e] & !dead[e]]; hold[j[neg]] <- H                     # a negative EGD starts the hold
    }
    # step 2: rest of the workup, path by path
    rest <- which(!found & !dead); pb_now <- pb[j[rest]]               # history at the start of this round
    for (b in 1:3) {
      rb <- rest[pb_now == b]; if (!length(rb)) next
      d <- PT[[b]]; isca <- type[rb] == "o"
      leaf <- integer(length(rb))
      if (any(isca)) leaf[isca] <- sample.int(nrow(d), sum(isca), TRUE, d$w_ca)
      if (any(!isca)) leaf[!isca] <- sample.int(nrow(d), sum(!isca), TRUE, d$w_free)
      L <- d[leaf, ]; dfp_v[rb] <- L$dfp
      for (kk in c("ent", "pet", "eus", "samp", "tb", "lap", "incid")) O[[kk]][j[rb]] <- O[[kk]][j[rb]] + L[[paste0("n_", kk)]]
      O$ww[j[rb]] <- O$ww[j[rb]] + L$n_ww
      nproc[rb] <- nproc[rb] + L$n_ent + L$n_pet + L$n_eus + L$n_incid + L$n_lap + L$n_ww
      days[rb] <- days[rb] + L$days
      sae <- runif(length(rb)) < L$sae; sdth <- sae & runif(length(rb)) < p$p_death_sae
      O$sae[j[rb]] <- O$sae[j[rb]] + sae
      cst[rb] <- cst[rb] + L$cost + sae * p$c_sae; qp[rb] <- qp[rb] + (sae & !sdth) * p$d_sae + L$rad * p$loss_rad_ca
      dp[rb] <- dp[rb] + L$rad * p$loss_rad_ca
      dead[rb] <- dead[rb] | sdth
      fo <- isca & runif(length(rb)) < L$p_found_ca & !sdth                  # another cancer found
      if (any(fo)) {
        pv <- (1 - pls) * V[[t]]["lead_early", ] + pls * V[[t]]["lead_late", ] - V[[t]]["pre", ]
        jj <- j[rb[fo]]; O$found_other[jj] <- O$found_other[jj] + 1
        O$cost[jj] <- O$cost[jj] + dfc * (pv[["cost"]] + (1 - pls) * p$c_tx_first_early + pls * p$c_tx_first_late)
        O$qaly[jj] <- O$qaly[jj] + dfq * pv[["qaly"]]; O$daly[jj] <- O$daly[jj] + dfd * pv[["daly"]]
      }
      fu <- type[rb] == "u" & L$n_pet > 0 & tau >= iv_until[j[rb]] & !sdth & runif(length(rb)) < (p$p_pet_ugi %||% 0)
      if (any(fu)) {                                                        # seen on PET -> targeted EGD
        found[rb[fu]] <- TRUE; O$egd[j[rb[fu]]] <- O$egd[j[rb[fu]]] + 1; cst[rb[fu]] <- cst[rb[fu]] + p$c_egd
        nproc[rb[fu]] <- nproc[rb[fu]] + 1
      }
      keep <- !(type[rb] == "u" & found[rb])
      pb[j[rb[keep]]] <- as.integer(L$flag[keep])
    }
    # valuation of the workup: false positives lose the whole-pathway disutility of where
    # their workup ended (Matza 2024; fp_tto = 1) or d_proc per procedure; people with
    # cancer lose d_proc per procedure
    fp <- type == "f"
    path_loss <- ifelse(fp & (p$fp_tto %||% 0) == 1, dfp_v, nproc * p$d_proc)
    path_daly <- ifelse(fp, days / 365.25 * (p$dw_anx %||% 0), nproc * p$d_proc)
    O$days[j] <- O$days[j] + days
    worry <- fp & !dead & runif(nj) < (p$p_worry %||% 0); O$worry[j] <- O$worry[j] + worry
    ww <- worry * (p$dur_worry %||% 0) * (p$dw_anx %||% 0)
    O$cost[j] <- O$cost[j] + dfc * cst
    O$qaly[j] <- O$qaly[j] - dfq * (qp + path_loss + ww)
    O$daly[j] <- O$daly[j] + dfd * (dp + (qp - dp) + path_daly + ww)
    O$cost_fp[j[fp]] <- O$cost_fp[j[fp]] + dfc * cst[fp]
    O$qaly_fp[j[fp]] <- O$qaly_fp[j[fp]] + dfq * (qp[fp] + path_loss[fp] + ww[fp])
    O$daly_fp[j[fp]] <- O$daly_fp[j[fp]] + dfd * (qp[fp] + path_daly[fp] + ww[fp])
    k_fp[j[fp]] <- k_fp[j[fp]] + 1L
    T_pd[j[dead]] <- tau
    fj <- j[found & !dead]; dx_t[fj] <- tau; dx_late[fj] <- runif(length(fj)) < pls
  }
  list(O = O, dx_t = dx_t, dx_late = dx_late, T_pd = T_pd, k_fp = k_fp, pm = pm)
}

# ---- Per-person timeline and outcomes -------------------------------------------
ms_outcomes <- function(nh, sc, prep, settings) {
  p <- prep$p; curve <- prep$curve; nT <- nrow(curve); N <- nrow(nh); O <- sc$O
  r_e <- -log(1 - p$m_early); r_l <- -log(1 - p$m_late)
  det <- is.finite(sc$dx_t)
  T_d0 <- pmin(nh$T_oc, sc$T_pd)                                   # death before presentation (other cause or procedure)
  clin <- !det & nh$T_clin < pmin(T_d0, nT)
  pres <- nh$T_clin < pmin(T_d0, nT)                               # the cancer would surface (screen-detected or not)
  st_late <- ifelse(det, sc$dx_late, nh$late)                      # stage from presentation on
  T_cd <- ifelse(pres, nh$T_clin - log(nh$u_cd) / ifelse(st_late, r_l, r_e), Inf)
  T_death <- pmin(T_d0, ifelse(pres, pmin(T_cd, nh$T_oc_post), Inf), nT)
  ca_death <- pres & T_cd <= pmin(nh$T_oc_post, nT) & T_cd < sc$T_pd
  proc_death <- is.finite(sc$T_pd) & sc$T_pd <= T_death
  D <- ifelse(det, sc$dx_t, ifelse(clin, nh$T_clin, Inf))          # diagnosis time
  lead_late <- det & sc$dx_late                                    # stage carried while in lead time
  ov <- function(a, b, t) pmax(0, pmin(b, t) - pmax(a, t - 1))    # time in [a, b) during year t
  ly <- qaly <- cost <- daly <- numeric(N)
  dq <- settings$disc_qaly; dc <- settings$disc_cost; dd <- settings$disc_daly
  for (t in seq_len(nT)) {
    alive <- ov(0, T_death, t)
    lead <- ov(ifelse(det, sc$dx_t, Inf), pmin(nh$T_clin, T_death), t)
    after <- ov(ifelse(pres, nh$T_clin, Inf), T_death, t)
    early <- ifelse(det & !lead_late, lead, 0) + ifelse(!st_late, after, 0)
    late <- ifelse(lead_late, lead, 0) + ifelse(st_late, after, 0)
    ly <- ly + alive
    qaly <- qaly + (1 + dq)^-(t - 1) * curve$utility_well[t] * (alive - early - late + p$u_mult_early * early + p$u_mult_late * late)
    cost <- cost + (1 + dc)^-(t - 1) * (p$c_fu_early * early + p$c_fu_late * late)
    daly <- daly + (1 + dd)^-(t - 1) * (p$dw_early * early + p$dw_late * late)
  }
  yr <- function(x) pmin(nT, floor(x) + 1)
  # symptomatic diagnosis (workup incl. an EGD that can perforate, first-year treatment)
  harm <- p$p_perf_egd * (1 - p$p_death_perf) * p$d_sae
  tc <- yr(nh$T_clin[clin])
  cost[clin] <- cost[clin] + (1 + dc)^-(tc - 1) * (p$c_symp_dx + ifelse(nh$late[clin], p$c_tx_first_late, p$c_tx_first_early) + p$p_perf_egd * p$c_perf)
  qaly[clin] <- qaly[clin] - (1 + dq)^-(tc - 1) * harm; daly[clin] <- daly[clin] + (1 + dd)^-(tc - 1) * harm
  # screen detection: first-year treatment at detection
  td <- yr(sc$dx_t[det])
  cost[det] <- cost[det] + (1 + dc)^-(td - 1) * ifelse(sc$dx_late[det], p$c_tx_first_late, p$c_tx_first_early)
  # deaths: end-of-life cost and years of life lost
  tk <- yr(T_cd[ca_death]); cost[ca_death] <- cost[ca_death] + (1 + dc)^-(tk - 1) * p$c_eol
  daly[ca_death] <- daly[ca_death] + (1 + dd)^-(tk - 1) * curve$daly_ref_le[tk]
  tp <- yr(sc$T_pd[proc_death]); daly[proc_death] <- daly[proc_death] + (1 + dd)^-(tp - 1) * curve$daly_ref_le[tp]
  data.frame(cancer = nh$T_on < pmin(T_death, nT), screen_detected = det, clin_dx = clin,
             overdiagnosed = det & !pres, stage_late = ifelse(det | clin, st_late, NA),
             onset_age = ifelse(is.finite(nh$T_on), nh$T_on, NA), dx_time = ifelse(is.finite(D), D, NA),
             ca_death = ca_death, proc_death = proc_death, death_time = T_death,
             ly = ly, qaly = qaly + O$qaly, cost = cost + O$cost, daly = daly + O$daly,
             tests = O$tests, positives = O$positives, false_pos = O$fp, fp_other_cancer = O$fp_other,
             egd = O$egd, egd_skipped = O$egd_skipped, ent = O$ent, pet = O$pet, eus = O$eus, samp = O$samp,
             tumor_board = O$tb, laparoscopy = O$lap, watchful_waiting = O$ww, incidental = O$incid,
             perforation = O$perf, sae = O$sae, days_waiting = O$days, hidden = O$hidden, worry = O$worry,
             found_other = O$found_other, cost_fp = O$cost_fp, qaly_fp = O$qaly_fp, daly_fp = O$daly_fp)
}

# One cohort row, all arms, the same simulated people (in chunks to bound memory).
ms_row <- function(row, data, p, settings, arms, n, seed, chunk = 2e5) {
  if (miscan_on(p)) stop("the MISCAN-ESO scenario (miscan_eac = 1) runs in the cohort model only; ",
                         "the microsimulation uses the base-case natural history")
  prep <- prepare_row(row, data, p, settings); wk <- workup_summary(prep$p)
  starts <- seq(1, n, by = chunk); out <- list(); inc <- NULL
  for (ci in seq_along(starts)) {
    m <- min(chunk, n - starts[ci] + 1)
    nh <- ms_natural_history(prep, m, seed + 7919L * ci, settings, row)
    for (a in seq_along(arms)) {
      sc <- ms_screen(nh, arms[a], prep, row, settings, wk, seed + 7919L * ci + 1000L * a)
      o <- ms_outcomes(nh, sc, prep, settings)
      o$person <- starts[ci] - 1 + seq_len(m); o$arm <- arms[a]; o$risk_group <- row$risk_group
      out[[length(out) + 1]] <- data.table::as.data.table(o)
      if (arms[a] == "Usual care") {                 # incidence check vs the SEER target
        nT <- nrow(prep$curve); Td <- o$death_time
        ev <- vapply(seq_len(nT), function(t) sum(o$clin_dx & o$dx_time >= t - 1 & o$dx_time < t), 0)
        al <- vapply(seq_len(nT), function(t) sum(Td > t - 1), 0)
        inc <- if (is.null(inc)) cbind(ev, al) else inc + cbind(ev, al)
      }
    }
    rm(nh); gc(FALSE)
  }
  res <- data.table::rbindlist(out)
  attr(res, "incidence") <- data.frame(age = row$start_age + seq_len(nrow(inc)) - 1, risk_group = row$risk_group,
    events = inc[, 1], alive = inc[, 2], target = target_ugi(prep$curve, row, settings))
  res
}
