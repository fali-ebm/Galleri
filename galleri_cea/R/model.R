# =============================================================
# PROJECT CONTRACT - Galleri upper GI screening & workup
# (engine in galleri_engine.R, figures in galleri_plots.R)
#
#   analyses$galleri_positive  people already Galleri +ve (esophagus/stomach
#                              signal): no workup vs EGD only vs full workup
#   analyses$screening         screening population: usual care vs Kunzmann
#                              score >= 8 -> EGD vs Galleri -> EGD only vs
#                              Galleri -> full workup
# Inputs: inputs/parameters.csv, baseline_risks.csv (competing risks by
# risk_group x sex x age), cohort.csv, counts.csv (observed workup counts).
# =============================================================

project_title <- "Galleri upper GI screening & workup"

settings <- ceakit::default_settings(
  wtp = 100000, n_sim = 1000, seed = 2026, owsa_pct = 0.25, evppi = TRUE,
  evppi_groups = list(
    `Galleri accuracy` = c("sens_g_ugi", "sens_g_other", "fpr_g_es"),
    `Natural history` = c("sojourn_years", "m_early", "m_late", "p_late_clin", "p_late_screen"),
    `Utilities & harms` = c("u_mult_early", "u_mult_late", "d_proc", "d_sae")),
  max_age = 100, disc_cost = 0.03, disc_qaly = 0.03, disc_daly = 0)

load_data <- function(dir) {
  rd <- function(stem) ceakit::read_table(list.files(dir, paste0("^", stem, "\\.(csv|xlsx|xls)$"), full.names = TRUE)[1])
  B <- rd("baseline_risks"); C <- rd("cohort")
  B$sex <- tolower(trimws(B$sex)); C$sex <- tolower(trimws(C$sex)); C$population <- tolower(trimws(C$population))
  list(baseline_risks = B, cohort = C)
}

validate_data <- function(data, settings) {
  e <- character(0); w <- character(0); B <- data$baseline_risks; C <- data$cohort
  need_b <- c("risk_group", "sex", "age_from", "age_to", "mort_other", "inc_ugi", "prev_ugi", "prev_other",
              "utility_well", "daly_ref_le")
  need_c <- c("population", "risk_group", "sex", "start_age", "share", "kunzmann_scoped")
  if (length(m <- setdiff(need_b, names(B)))) return(list(errors = paste("baseline_risks: missing", paste(m, collapse = ", ")), warnings = w))
  if (length(m <- setdiff(need_c, names(C)))) return(list(errors = paste("cohort: missing", paste(m, collapse = ", ")), warnings = w))
  for (cc in c("mort_other", "inc_ugi", "prev_ugi", "prev_other", "utility_well"))
    if (any(is.na(B[[cc]]) | B[[cc]] < 0 | B[[cc]] > 1)) e <- c(e, sprintf("baseline_risks: %s must be 0-1", cc))
  if (any(B$prev_ugi + B$prev_other > 1, na.rm = TRUE)) e <- c(e, "baseline_risks: prev_ugi + prev_other > 1")
  for (pp in c("screening", "galleri_positive")) {
    s <- sum(C$share[C$population == pp])
    if (abs(s - 1) > 0.001) e <- c(e, sprintf("cohort: shares for '%s' sum to %.4f (must be 1)", pp, s))
  }
  if (!all(C$kunzmann_scoped %in% c(0, 1))) e <- c(e, "cohort: kunzmann_scoped must be 0 or 1")
  if (any(C$start_age >= settings$max_age)) e <- c(e, "cohort: start_age must be below max_age")
  if (!length(e)) for (i in seq_len(nrow(C))) tryCatch(baseline_for(data, C[i, ], settings),
    error = function(err) e <<- c(e, sprintf("cohort row %d: %s", i, conditionMessage(err))))
  if (any(toupper(c(B$status, C$status)) == "DUMMY", na.rm = TRUE)) w <- "baseline_risks / cohort still contain DUMMY rows"
  list(errors = e, warnings = w)
}

METRICS <- tibble::tribble(
  ~key, ~section, ~label, ~type, ~digits,
  "total_cancers", "1. Baseline burden", "Prevalent cancers at screening (same in every arm)", "count", 1,
  "cancer_free", "1. Baseline burden", "Cancer-free people", "count", 0,
  "kunz_scores", "2. Screening & workup burden", "Kunzmann scores", "count", 0,
  "galleri_tests", "2. Screening & workup burden", "Galleri tests", "count", 0,
  "investigated", "2. Screening & workup burden", "People investigated (EGD or more)", "count", 0,
  "egd", "2. Screening & workup burden", "EGDs", "count", 0,
  "ent", "2. Screening & workup burden", "ENT evaluations", "count", 0,
  "pet_ct", "2. Screening & workup burden", "PET CT scans", "count", 0,
  "eus", "2. Screening & workup burden", "EUS procedures", "count", 0,
  "eus_sampling", "2. Screening & workup burden", "EUS samplings", "count", 0,
  "tumor_board", "2. Screening & workup burden", "Tumor board discussions", "count", 0,
  "invasive_total", "2. Screening & workup burden", "Invasive procedures (EGD/EUS/sampling)", "count", 0,
  "mean_days", "2. Screening & workup burden", "Mean days to resolution if investigated", "ratio", 0,
  "cost_tests", "2. Screening & workup burden", "Screening test cost per person ($)", "person", 2,
  "cost_dx", "2. Screening & workup burden", "Diagnostic workup cost per person ($)", "person", 2,
  "cost_tx", "2. Screening & workup burden", "Lifetime treatment & care cost per person ($)", "person", 2,
  "cost_total", "2. Screening & workup burden", "Total lifetime cost per person ($)", "person", 2,
  "screen_detected", "3. Benefit", "Cancers found by screening/workup", "count", 1,
  "found_upper_gi", "3. Benefit", "  - upper GI (esophagus/stomach)", "count", 1,
  "found_other", "3. Benefit", "  - other sites", "count", 1,
  "clin_dx_lifetime", "3. Benefit", "Cancers diagnosed with symptoms over lifetime", "count", 1,
  "late_lifetime", "3. Benefit", "Late-stage diagnoses over lifetime", "count", 1,
  "late_averted", "3. Benefit", "Late-stage diagnoses averted vs reference", "count", 1,
  "ca_deaths", "3. Benefit", "Cancer deaths over lifetime", "count", 1,
  "ca_deaths_averted", "3. Benefit", "Cancer deaths averted vs reference", "count", 1,
  "oc_deaths", "3. Benefit", "Other-cause deaths (competing risk)", "count", 0,
  "ly_total", "3. Benefit", "Life years (undiscounted)", "count", 0,
  "ly_gained", "3. Benefit", "Life years gained vs reference", "count", 1,
  "free_investigated", "4. Harm", "Cancer-free people investigated (false positives)", "count", 0,
  "invasive_free", "4. Harm", "Invasive procedures in cancer-free people", "count", 0,
  "sae_total", "4. Harm", "Serious procedure complications", "count", 2,
  "sae_free", "4. Harm", "  - in cancer-free people", "count", 2,
  "proc_deaths", "4. Harm", "Deaths from procedure complications", "count", 3,
  "dose_sv", "4. Harm", "Collective PET CT dose (Sv)", "count", 2,
  "radiation_cancers", "4. Harm", "Radiation-induced cancers (lifetime)", "count", 3,
  "missed_total", "4. Harm", "Prevalent cancers missed at screening", "count", 1,
  "missed_investigated", "4. Harm", "  - missed despite investigation (false reassurance)", "count", 1,
  "qaly_total", "5. QALYs & DALYs", "Lifetime QALYs (discounted)", "count", 0,
  "qaly_harms", "5. QALYs & DALYs", "QALYs lost to procedure harms", "count", 1,
  "qaly_net", "5. QALYs & DALYs", "Net QALYs gained vs reference", "count", 1,
  "daly_total", "5. QALYs & DALYs", "Cancer DALYs (YLL + YLD + harms)", "count", 0,
  "daly_harms", "5. QALYs & DALYs", "  - DALYs from procedure harms", "count", 1,
  "daly_averted", "5. QALYs & DALYs", "Net DALYs averted vs reference", "count", 1,
  "cost_per_daly", "6. Efficiency", "Cost per DALY averted vs reference ($)", "ratio", 0,
  "nn_per_found", "6. Efficiency", "People tested/worked up per cancer found", "ratio", 1,
  "nn_per_ca_death_averted", "6. Efficiency", "People tested/worked up per cancer death averted", "ratio", 0,
  "investigated_per_found", "6. Efficiency", "People investigated per cancer found", "ratio", 1,
  "fp_per_found", "6. Efficiency", "False positives investigated per cancer found", "ratio", 1,
  "invasive_per_found", "6. Efficiency", "Invasive procedures per cancer found", "ratio", 1)

ALL_PARAMS <- c("sens_g_ugi", "sens_g_other", "fpr_g_es", "sens_egd", "sens_ds_other", "p_egd_pos",
  "p_ent_route", "p_ent_pos", "p_pet_pb", "p_pet_incid", "p_eus_pos", "p_eus_neg", "p_miss",
  "p_incid_cancer", "p_late_screen", "p_late_clin", "sojourn_years", "m_early", "m_late", "sae_egd",
  "sae_ent", "sae_eus", "sae_eus_samp", "p_death_sae", "msv_pet", "risk_per_sv", "loss_rad_ca", "t_egd",
  "t_ent", "t_pet", "t_eus", "t_tb", "c_galleri", "c_kunz", "c_egd", "c_ent", "c_pet", "c_eus",
  "c_eus_samp", "c_tb", "c_symp_dx", "c_tx_first_early", "c_tx_first_late", "c_fu_early", "c_fu_late",
  "c_eol", "c_sae", "u_mult_early", "u_mult_late", "d_proc", "d_sae", "dw_early", "dw_late")

keep_ordered <- function(q) {   # sampled values stay logically ordered
  q$u_mult_late <- min(q$u_mult_late, q$u_mult_early); q
}

KEY <- c("screen_detected", "late_averted", "ca_deaths_averted", "free_investigated",
         "invasive_free", "sae_total", "qaly_net", "daly_averted")

make_analysis <- function(name, title, arms, owsa, cohort_n, unit, required, plots)
  ceakit::cea_analysis(
    name = name, title = title, arms = arms,
    run = function(p, data, settings) run_population(p, data, settings, name)$total,
    cost = "cost_total", effect = "qaly_total", effect_label = "QALY",
    metrics = METRICS, derive = derive_galleri, required = required, cohort_n = cohort_n, unit = unit,
    owsa = owsa, key_metrics = KEY, constraints = keep_ordered, plots = plots)

analyses <- list(
  galleri_positive = make_analysis(
    "galleri_positive", "Galleri +ve workup", ARMS_GP, c("Full workup", "No workup"), 1000,
    "per 1,000 Galleri +ve",
    setdiff(ALL_PARAMS, c("sens_g_ugi", "sens_g_other", "fpr_g_es", "sens_egd", "sens_ds_other", "c_galleri", "c_kunz")),
    list(workup_tree = plot_workup_tree, benefit_harm = plot_balance, icon_array = plot_icons,
         markov_trace = plot_trace, cancer_states = plot_cancer_states, subgroups = plot_subgroups)),
  screening = make_analysis(
    "screening", "Screening comparison", ARMS_SC, c("Galleri -> full workup", "Kunzmann -> EGD"), 100000,
    "per 100,000 screened",
    setdiff(ALL_PARAMS, c("p_egd_pos", "p_ent_pos", "p_pet_incid", "p_miss", "p_incid_cancer")),
    list(strategy_tree = plot_strategy_tree, benefit_harm = plot_balance, markov_trace = plot_trace,
         cancer_states = plot_cancer_states, subgroups = plot_subgroups)))
