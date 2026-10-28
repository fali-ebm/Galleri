# =============================================================
# Galleri model engine (project-specific). Uses ceakit engines:
#   ceakit::tree_paths()     - Galleri +ve workup decision tree
#   ceakit::markov_values()  - lifetime Markov with competing risks
#   ceakit::lookup_by_age()  - baseline risks by risk group x sex x age
# =============================================================

# lead_early / lead_late: screen-detected cancers during their lead time (the
# time until they would have presented with symptoms). No cancer deaths occur
# in lead time, so screening is credited only for stage shift, not for the
# lead time itself (stage-specific survival is measured from clinical diagnosis).
STATES <- c("well", "pre", "lead_early", "lead_late", "early", "late", "dead_cancer", "dead_other")
DEAD <- c("dead_cancer", "dead_other")
CHANNELS <- c("cost", "qaly", "daly", "ly", "clin_dx", "clin_late", "ca_death", "oc_death")
ARMS_GP <- c("No workup", "EGD only", "Full workup")
ARMS_SC <- c("No screening (usual care)", "Kunzmann -> EGD", "Galleri -> EGD only", "Galleri -> full workup")

# ---- Markov model ---------------------------------------------------------
baseline_for <- function(data, row, settings)
  ceakit::lookup_by_age(data$baseline_risks, list(risk_group = row$risk_group, sex = row$sex),
                        row$start_age:(settings$max_age - 1))

transition_matrix <- function(mo, inc, p) {
  pp <- 1 - exp(-1 / p$sojourn_years); pl <- p$p_late_clin; s <- 1 - mo
  n <- length(STATES)
  P <- matrix(0, n, n, dimnames = list(STATES, STATES))
  P[, "dead_other"] <- mo
  P["well", c("well", "pre")] <- s * c(1 - inc, inc)
  P["pre", c("pre", "early", "late")] <- s * c(1 - pp, pp * (1 - pl), pp * pl)
  P["lead_early", c("lead_early", "early")] <- s * c(1 - pp, pp)
  P["lead_late", c("lead_late", "late")] <- s * c(1 - pp, pp)
  P["early", c("early", "dead_cancer")] <- s * c(1 - p$m_early, p$m_early)
  P["late", c("late", "dead_cancer")] <- s * c(1 - p$m_late, p$m_late)
  P[DEAD, ] <- 0; P["dead_cancer", "dead_cancer"] <- 1; P["dead_other", "dead_other"] <- 1
  P
}

markov_parts <- function(curve, p) {
  transition <- function(t) transition_matrix(curve$mort_other[t], curve$inc_ugi[t], p)
  rewards <- function(t, P) {
    alive <- setNames(ceakit::alive_fraction(P, DEAD), STATES)
    live <- !STATES %in% DEAD
    R <- matrix(0, length(STATES), length(CHANNELS), dimnames = list(STATES, CHANNELS))
    u <- c(well = 1, pre = 1, lead_early = p$u_mult_early, lead_late = p$u_mult_late,
           early = p$u_mult_early, late = p$u_mult_late, dead_cancer = 0, dead_other = 0)
    dw <- c(well = 0, pre = 0, lead_early = p$dw_early, lead_late = p$dw_late,
            early = p$dw_early, late = p$dw_late, dead_cancer = 0, dead_other = 0)
    R[, "ly"] <- alive
    R[, "qaly"] <- alive * curve$utility_well[t] * u[STATES]
    R["lead_early", "cost"] <- alive[["lead_early"]] * p$c_fu_early
    R["lead_late", "cost"]  <- alive[["lead_late"]] * p$c_fu_late
    R["early", "cost"] <- alive[["early"]] * p$c_fu_early + P["early", "dead_cancer"] * p$c_eol
    R["late", "cost"]  <- alive[["late"]] * p$c_fu_late  + P["late", "dead_cancer"]  * p$c_eol
    R["pre", "cost"]   <- P["pre", "early"] * (p$c_symp_dx + p$c_tx_first_early) +
                          P["pre", "late"]  * (p$c_symp_dx + p$c_tx_first_late)
    R[, "daly"] <- alive * dw[STATES] +                                       # YLD
                   P[, "dead_cancer"] * curve$daly_ref_le[t] * live           # YLL, new deaths only
    R["pre", "clin_dx"] <- P["pre", "early"] + P["pre", "late"]
    R["pre", "clin_late"] <- P["pre", "late"]
    R[, "ca_death"] <- P[, "dead_cancer"] * live
    R[, "oc_death"] <- P[, "dead_other"] * live
    R
  }
  list(n = nrow(curve), transition = transition, rewards = rewards)
}

galleri_markov_values <- function(curve, p, settings) {
  mk <- markov_parts(curve, p)
  ceakit::markov_values(mk$n, mk$transition, mk$rewards,
                        disc = c(cost = settings$disc_cost, qaly = settings$disc_qaly, daly = settings$disc_daly))
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
  n_tb   = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1))

workup_paths <- function(p) {
  d <- ceakit::tree_paths(workup_edges(p), WORKUP_NODES)
  found <- c(egd_pos = 1, ent_pos = 1, eus_pos = 1, incid = p$p_incid_cancer)
  site <- c(egd_pos = "Upper GI", ent_pos = "Head & neck", eus_pos = "Pancreatobiliary", incid = "Other (incidental)")
  d$found <- ifelse(d$leaf %in% names(found), found[d$leaf], 0)
  d$missed <- ifelse(d$leaf %in% c("ent_neg", "eus_neg", "no_eus", "pet_neg"), p$p_miss, 0)
  d$site <- unname(site[d$leaf])
  add_path_details(d, p)
}

egd_only_paths <- function(p, full = workup_paths(p)) {
  rest <- full[full$leaf != "egd_pos", ]
  ca_neg <- sum(rest$prob * (rest$found + rest$missed)) / sum(rest$prob)
  tibble::tibble(leaf = c("egd_pos", "egd_neg"), prob = c(p$p_egd_pos, 1 - p$p_egd_pos),
                 path = c("EGD +ve", "EGD -ve > stop"), n_egd = 1, n_ent = 0, n_pet = 0, n_eus = 0,
                 n_samp = 0, n_tb = 0, found = c(1, 0), missed = c(0, ca_neg), site = c("Upper GI", NA)) |>
    add_path_details(p)
}

add_path_details <- function(d, p) {
  d$free <- 1 - d$found - d$missed
  d$n_invasive <- d$n_egd + d$n_eus + d$n_samp
  d$sae <- d$n_egd * p$sae_egd + d$n_ent * p$sae_ent + d$n_eus * p$sae_eus + d$n_samp * p$sae_eus_samp
  d$rad_ca <- d$n_pet * p$msv_pet / 1000 * p$risk_per_sv
  d$days <- d$n_egd * p$t_egd + d$n_ent * p$t_ent + d$n_pet * p$t_pet + d$n_eus * p$t_eus + d$n_tb * p$t_tb
  d$cost_dx <- d$n_egd * p$c_egd + d$n_ent * p$c_ent + d$n_pet * p$c_pet + d$n_eus * p$c_eus +
               d$n_samp * p$c_eus_samp + d$n_tb * p$c_tb
  d
}

# ---- Ledgers: events at time 0, per person ---------------------------------
ledger <- function(total) list(
  total = total, found = 0, found_ugi = 0, missed = total, missed_inv = 0, investigated = 0,
  free_inv = 0, kunz = 0, galleri = 0, n_egd = 0, n_ent = 0, n_pet = 0, n_eus = 0, n_samp = 0,
  n_tb = 0, n_invasive = 0, invasive_free = 0, sae = 0, sae_free = 0, rad_ca = 0,
  cost_tests = 0, cost_dx = 0, days = 0)

ledger_from_paths <- function(d, total, scale = 1, galleri = 0, test_cost = 0) {
  s <- function(x) scale * sum(d$prob * x)
  l <- ledger(total)
  l$found <- s(d$found); l$found_ugi <- s(d$found * (d$site %in% "Upper GI"))
  l$missed <- total - l$found; l$missed_inv <- s(d$missed)
  l$investigated <- scale; l$free_inv <- s(d$free)
  for (k in c("n_egd", "n_ent", "n_pet", "n_eus", "n_samp", "n_tb", "n_invasive", "sae", "rad_ca",
              "cost_dx", "days")) l[[k]] <- s(d[[k]])
  l$invasive_free <- s(d$free * d$n_invasive); l$sae_free <- s(d$free * d$sae)
  l$galleri <- galleri; l$cost_tests <- test_cost
  l
}

arms_galleri_positive <- function(p, row = NULL, curve0 = NULL) {
  full <- workup_paths(p); egd <- egd_only_paths(p, full)
  total <- sum(full$prob * (full$found + full$missed))
  stats::setNames(list(ledger(total), ledger_from_paths(egd, total), ledger_from_paths(full, total)), ARMS_GP)
}

arms_screening <- function(p, row, curve0) {
  pu <- curve0$prev_ugi; po <- curve0$prev_other; total <- pu + po; free <- 1 - total
  a_u <- pu * p$sens_g_ugi; a_o <- po * p$sens_g_other; a_f <- free * p$fpr_g_es; gp <- a_u + a_o + a_f
  galleri <- function(full) {
    l <- ledger(total); l$galleri <- 1; l$cost_tests <- p$c_galleri
    l$investigated <- gp; l$free_inv <- a_f; l$n_egd <- gp
    found_u <- a_u * p$sens_egd
    if (!full) {
      l$found <- l$found_ugi <- found_u; l$missed_inv <- a_u - found_u + a_o
      l$n_invasive <- gp; l$invasive_free <- a_f
      l$sae <- gp * p$sae_egd; l$sae_free <- a_f * p$sae_egd
      l$cost_dx <- gp * p$c_egd; l$days <- gp * p$t_egd
    } else {
      n2 <- gp - found_u; found_o <- a_o * p$sens_ds_other
      ent <- n2 * p$p_ent_route; pet <- n2 * (1 - p$p_ent_route)
      eus <- pet * p$p_pet_pb; samp <- eus * (p$p_eus_pos + p$p_eus_neg)
      tb <- ent * (1 - p$p_ent_pos) + pet * (1 - p$p_pet_pb * p$p_eus_pos)
      f2 <- if (n2 > 0) a_f / n2 else 0
      l$found <- found_u + found_o; l$found_ugi <- found_u
      l$missed_inv <- a_u - found_u + a_o - found_o
      l$n_ent <- ent; l$n_pet <- pet; l$n_eus <- eus; l$n_samp <- samp; l$n_tb <- tb
      l$n_invasive <- gp + eus + samp; l$invasive_free <- a_f + f2 * (eus + samp)
      ds <- ent * p$sae_ent + eus * p$sae_eus + samp * p$sae_eus_samp
      l$sae <- gp * p$sae_egd + ds; l$sae_free <- a_f * p$sae_egd + f2 * ds
      l$rad_ca <- pet * p$msv_pet / 1000 * p$risk_per_sv
      l$cost_dx <- gp * p$c_egd + ent * p$c_ent + pet * p$c_pet + eus * p$c_eus + samp * p$c_eus_samp + tb * p$c_tb
      l$days <- gp * p$t_egd + ent * p$t_ent + pet * p$t_pet + eus * p$t_eus + tb * p$t_tb
    }
    l$missed <- total - l$found
    l
  }
  k <- ledger(total); k$kunz <- 1; k$cost_tests <- p$c_kunz
  if (isTRUE(row$kunzmann_scoped == 1)) {
    k$found <- k$found_ugi <- pu * p$sens_egd; k$missed <- total - k$found
    k$missed_inv <- total - k$found; k$investigated <- 1; k$free_inv <- free
    k$n_egd <- k$n_invasive <- 1; k$invasive_free <- free
    k$sae <- p$sae_egd; k$sae_free <- free * p$sae_egd; k$cost_dx <- p$c_egd; k$days <- p$t_egd
  }
  stats::setNames(list(ledger(total), k, galleri(FALSE), galleri(TRUE)), ARMS_SC)
}

# ---- Outcomes per cohort row and per population ----------------------------
ADDITIVE <- c("total_cancers", "cancer_free", "kunz_scores", "galleri_tests", "investigated",
  "egd", "ent", "pet_ct", "eus", "eus_sampling", "tumor_board", "invasive_total", "days_x_inv",
  "cost_tests", "cost_dx", "cost_total", "screen_detected", "found_upper_gi", "found_other",
  "late_lifetime", "clin_dx_lifetime", "ca_deaths", "oc_deaths", "ly_total", "free_investigated",
  "invasive_free", "sae_total", "sae_free", "proc_deaths", "dose_sv", "radiation_cancers",
  "missed_total", "missed_investigated", "qaly_total", "qaly_harms", "daly_total", "daly_harms")

initial_state <- function(l, p) {
  fl <- l$found * p$p_late_screen; pd <- l$sae * p$p_death_sae
  c(well = max(0, 1 - l$total - pd), pre = l$missed, lead_early = l$found - fl, lead_late = fl,
    early = 0, late = 0, dead_cancer = 0, dead_other = pd)
}

row_outcomes <- function(row, data, p, settings, population) {
  curve <- baseline_for(data, row, settings)
  V <- galleri_markov_values(curve, p, settings)
  arms <- if (population == "screening") arms_screening(p, row, curve[1, ]) else arms_galleri_positive(p)
  one <- function(l) {
    init <- initial_state(l, p); life <- drop(init %*% V)
    fl <- init[["lead_late"]]; pd <- init[["dead_other"]]
    n_proc <- l$n_egd + l$n_ent + l$n_pet + l$n_eus
    q_harm <- n_proc * p$d_proc + (l$sae - pd) * p$d_sae + l$rad_ca * p$loss_rad_ca
    d_harm <- q_harm + pd * curve$daly_ref_le[1]
    cost0 <- l$cost_tests + l$cost_dx + (l$found - fl) * p$c_tx_first_early + fl * p$c_tx_first_late + l$sae * p$c_sae
    c(total_cancers = l$total, cancer_free = 1 - l$total, kunz_scores = l$kunz, galleri_tests = l$galleri,
      investigated = l$investigated, egd = l$n_egd, ent = l$n_ent, pet_ct = l$n_pet, eus = l$n_eus,
      eus_sampling = l$n_samp, tumor_board = l$n_tb, invasive_total = l$n_invasive, days_x_inv = l$days,
      cost_tests = l$cost_tests, cost_dx = l$cost_dx, cost_total = cost0 + life[["cost"]],
      screen_detected = l$found, found_upper_gi = l$found_ugi, found_other = l$found - l$found_ugi,
      late_lifetime = fl + life[["clin_late"]], clin_dx_lifetime = life[["clin_dx"]],
      ca_deaths = life[["ca_death"]], oc_deaths = life[["oc_death"]], ly_total = life[["ly"]],
      free_investigated = l$free_inv, invasive_free = l$invasive_free, sae_total = l$sae,
      sae_free = l$sae_free, proc_deaths = pd, dose_sv = l$n_pet * p$msv_pet / 1000,
      radiation_cancers = l$rad_ca, missed_total = l$missed, missed_investigated = l$missed_inv,
      qaly_total = life[["qaly"]] - q_harm, qaly_harms = q_harm,
      daly_total = life[["daly"]] + d_harm, daly_harms = d_harm)[ADDITIVE]
  }
  sapply(arms, one)
}

cohort_rows <- function(data, population) data$cohort[data$cohort$population == population, ]

run_population <- function(p, data, settings, population) {
  C <- cohort_rows(data, population)
  mats <- lapply(seq_len(nrow(C)), function(i) row_outcomes(C[i, ], data, p, settings, population))
  list(cohort = C, rows = mats, total = Reduce(`+`, Map(`*`, mats, C$share)))
}

# Galleri-specific derived metrics (ceakit adds inc_cost, inc_effect, icer, nmb)
derive_galleri <- function(A, settings) {
  per <- function(x, y) ifelse(y > 1e-15, x / y, NA)
  r <- A[, 1]; not_ref <- seq_len(ncol(A)) != 1
  d_av <- r["daly_total"] - A["daly_total", ]; ca_av <- r["ca_deaths"] - A["ca_deaths", ]
  rbind(mean_days = per(A["days_x_inv", ], A["investigated", ]),
        cost_tx = A["cost_total", ] - A["cost_tests", ] - A["cost_dx", ],
        late_averted = r["late_lifetime"] - A["late_lifetime", ],
        ca_deaths_averted = ca_av, ly_gained = A["ly_total", ] - r["ly_total"],
        qaly_net = A["qaly_total", ] - r["qaly_total"], daly_averted = d_av,
        cost_per_daly = ifelse(not_ref, (A["cost_total", ] - r["cost_total"]) / d_av, NA),
        nn_per_found = per(1, A["screen_detected", ]),
        nn_per_ca_death_averted = ifelse(not_ref, per(1, ca_av), NA),
        investigated_per_found = per(A["investigated", ], A["screen_detected", ]),
        fp_per_found = per(A["free_investigated", ], A["screen_detected", ]),
        invasive_per_found = per(A["invasive_total", ], A["screen_detected", ]))
}

# Subgroup (cohort row) results, base case
subgroup_table <- function(p, data, settings, population, cohort_n) {
  pop <- run_population(p, data, settings, population)
  purrr::map_dfr(seq_along(pop$rows), function(i) {
    m <- rbind(pop$rows[[i]], derive_galleri(pop$rows[[i]], settings)); r <- pop$cohort[i, ]
    tibble::tibble(risk_group = r$risk_group, sex = r$sex, start_age = r$start_age, share = r$share,
                   arm = colnames(m), prevalent_cancers = m["total_cancers", ] * cohort_n,
                   cancers_found = m["screen_detected", ] * cohort_n,
                   false_positives = m["free_investigated", ] * cohort_n,
                   cancer_deaths_averted = m["ca_deaths_averted", ] * cohort_n,
                   net_qalys = m["qaly_net", ] * cohort_n, net_dalys_averted = m["daly_averted", ] * cohort_n,
                   incremental_cost_per_person = m["cost_total", ] - m["cost_total", 1])
  })
}

# Cohort-weighted Markov trace by arm
trace_table <- function(p, data, settings, population) {
  C <- cohort_rows(data, population)
  purrr::map_dfr(seq_len(nrow(C)), function(i) {
    row <- C[i, ]; curve <- baseline_for(data, row, settings); mk <- markov_parts(curve, p)
    arms <- if (population == "screening") arms_screening(p, row, curve[1, ]) else arms_galleri_positive(p)
    purrr::imap_dfr(arms, function(l, a) ceakit::markov_trace(mk$n, mk$transition, initial_state(l, p)) |>
      dplyr::mutate(arm = a, share = row$share))
  }) |>
    tidyr::pivot_longer(dplyr::all_of(STATES), names_to = "state", values_to = "prop") |>
    dplyr::mutate(state = sub("^lead_", "", state)) |>
    dplyr::group_by(arm, cycle, state) |> dplyr::summarise(prop = sum(prop * share), .groups = "drop") |>
    dplyr::mutate(state = factor(state, setdiff(STATES, c("lead_early", "lead_late")), c("Well", "Undiagnosed cancer", "Early-stage dx",
                                                   "Late-stage dx", "Cancer death", "Other-cause death")))
}
