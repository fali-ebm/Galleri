# =============================================================
# PROJECT CONTRACT - Galleri upper GI screening & workup
# (engine in galleri_engine.R, figures in galleri_plots.R)
#
# Annual screening (Galleri every year; EGD every year for people with a
# Kunzmann score >= 8) from the start age to settings$screen_stop_age, lifetime
# horizon. Natural history is calibrated to SEER incidence for each subgroup and
# drives every arm. Base case: 50-year-old non-Hispanic White men (cohort.csv).
#   analyses$galleri     usual care vs Galleri: benefits and harms per 100,000
#   analyses$comparison  usual care vs Kunzmann -> EGD vs Galleri
# Subgroup sensitivity grid: start age x sex x race/ethnicity (settings$subgroups).
# Inputs: inputs/parameters.csv, baseline_risks.csv (race/ethnicity x sex x
# age: life tables + SEER targets), cohort.csv, counts.csv (workup counts).
# =============================================================

project_title <- "Galleri upper GI screening & workup"

# Defaults; inputs/settings.csv (written by the workbook importer or by hand)
# overrides any of the scalar settings below when the project is loaded.
SETTINGS_DEFAULTS <- list(
  wtp = 100000, n_sim = 1000, seed = 2026, owsa_pct = 0.25, evppi = TRUE,
  evppi_groups = list(
    `Galleri accuracy` = c("sens_g_ugi", "sens_g_other", "fpr_g_es", "sens_g_rep_ratio"),
    `Natural history` = c("sojourn_years", "m_early", "m_late", "or_late_clin", "p_late_screen"),
    `Utilities & harms` = c("u_mult_early", "u_mult_late", "d_proc", "d_sae")),
  max_age = 100, disc_cost = 0.03, disc_qaly = 0.03, disc_daly = 0,
  screen_stop_age = 85, screen_interval = 1,          # annual screening to age 85
  # Kunzmann score >= 8 by sex (Kunzmann 2018: 59.7% of men, 3.5% of women flagged;
  # capture within sex derived assuming equal within-sex RR and 85% male cases)
  kunz_share_male = 0.597, kunz_share_female = 0.035, kunz_capture_male = 0.881, kunz_capture_female = 0.153,
  cost_year = 2025,                                    # costs in inputs are 2025 USD (data-raw/inflate_costs.R)
  subgroups = list(ages = seq(50, 90, 5), sexes = c("male", "female"),
                   races = c("all", "nhw", "nhb", "nhapi", "nhaian", "hisp")))
make_settings <- function(over = list()) do.call(ceakit::default_settings, utils::modifyList(SETTINGS_DEFAULTS, over))
settings <- make_settings()

# inputs/settings.csv: columns setting, value (other columns are documentation).
read_settings_file <- function(dir) {
  f <- file.path(dir, "settings.csv")
  if (!file.exists(f)) return(list())
  s <- utils::read.csv(f, colClasses = "character", check.names = FALSE)
  v <- stats::setNames(trimws(s$value), trimws(s$setting))
  lst <- function(k) trimws(strsplit(v[[k]], ",")[[1]])
  out <- list(); sub <- SETTINGS_DEFAULTS$subgroups
  for (k in names(v)) {
    if (k == "subgroup_ages") sub$ages <- as.numeric(lst(k))
    else if (k == "subgroup_sexes") sub$sexes <- tolower(lst(k))
    else if (k == "subgroup_races") sub$races <- tolower(lst(k))
    else if (k %in% names(SETTINGS_DEFAULTS)) out[[k]] <- as.numeric(v[[k]])
    else stop("inputs/settings.csv: unknown setting '", k, "'", call. = FALSE)
  }
  bad <- names(out)[vapply(out, function(x) length(x) != 1 || is.na(x), TRUE)]
  if (length(bad) || anyNA(sub$ages)) stop("inputs/settings.csv: non-numeric value for ",
                                          paste(c(bad, if (anyNA(sub$ages)) "subgroup_ages"), collapse = ", "), call. = FALSE)
  out$subgroups <- sub
  out
}

load_data <- function(dir) {
  rd <- function(stem) ceakit::read_table(list.files(dir, paste0("^", stem, "\\.(csv|xlsx|xls)$"), full.names = TRUE)[1])
  B <- rd("baseline_risks"); C <- rd("cohort")
  B$sex <- tolower(trimws(B$sex)); C$sex <- tolower(trimws(C$sex)); C$population <- tolower(trimws(C$population))
  B$race_eth <- tolower(trimws(B$race_eth)); C$race_eth <- tolower(trimws(C$race_eth))
  # settings.csv is read here because load_data is the hook that knows the inputs folder
  settings <<- make_settings(read_settings_file(dir))
  list(baseline_risks = B, cohort = C)
}

validate_data <- function(data, settings) {
  e <- character(0); w <- character(0); B <- data$baseline_risks; C <- data$cohort
  need_b <- c("race_eth", "sex", "age_from", "age_to", "mort_other", "inc_eso", "inc_stom", "inc_other",
              "eac_share", "late_eso", "late_stom", "utility_well", "daly_ref_le")
  need_c <- c("population", "race_eth", "risk_group", "sex", "start_age", "share", "kunzmann_scoped")
  if (length(m <- setdiff(need_b, names(B)))) return(list(errors = paste("baseline_risks: missing", paste(m, collapse = ", ")), warnings = w))
  if (length(m <- setdiff(need_c, names(C)))) return(list(errors = paste("cohort: missing", paste(m, collapse = ", ")), warnings = w))
  for (cc in c("mort_other", "inc_eso", "inc_stom", "inc_other", "eac_share", "late_eso", "late_stom", "utility_well"))
    if (any(is.na(B[[cc]]) | B[[cc]] < 0 | B[[cc]] > 1)) e <- c(e, sprintf("baseline_risks: %s must be 0-1", cc))
  miss <- setdiff(settings$subgroups$races, unique(B$race_eth))
  if (length(miss)) e <- c(e, paste("baseline_risks: no rows for subgroup race/ethnicity", paste(miss, collapse = ", ")))
  for (pp in POPULATION) {
    s <- sum(C$share[C$population == pp])
    if (abs(s - 1) > 0.001) e <- c(e, sprintf("cohort: shares for '%s' sum to %.4f (must be 1)", pp, s))
  }
  if (!all(C$kunzmann_scoped %in% c(0, 1))) e <- c(e, "cohort: kunzmann_scoped must be 0 or 1")
  for (i in which(C$risk_group %in% c("kunzmann_high", "kunzmann_low"))) {
    k <- settings[[paste0("kunz_share_", C$sex[i])]]
    want <- if (C$risk_group[i] == "kunzmann_high") k else 1 - k
    if (is.null(k) || abs(C$share[i] - want) > 1e-6)
      e <- c(e, sprintf("cohort row %d: %s share %.4f must equal settings kunz_share_%s (%s)", i, C$risk_group[i],
                        C$share[i], C$sex[i], format(want)))
  }
  if (any(C$start_age >= settings$max_age)) e <- c(e, "cohort: start_age must be below max_age")
  if (!length(e)) for (i in seq_len(nrow(C))) tryCatch(baseline_for(data, C[i, ], settings),
    error = function(err) e <<- c(e, sprintf("cohort row %d: %s", i, conditionMessage(err))))
  if (any(toupper(c(B$status, C$status)) == "DUMMY", na.rm = TRUE)) w <- "baseline_risks / cohort still contain DUMMY rows"
  list(errors = e, warnings = w)
}

METRICS <- tibble::tribble(
  ~key, ~section, ~label, ~type, ~digits,
  "total_cancers", "1. Baseline burden", "Preclinical upper GI cancers at the start age (calibrated)", "count", 1,
  "cancer_free", "1. Baseline burden", "People without preclinical upper GI cancer at the start", "count", 0,
  "kunz_scores", "2. Screening & workup burden", "Kunzmann scores", "count", 0,
  "galleri_tests", "2. Screening & workup burden", "Galleri tests (all screening rounds)", "count", 0,
  "screen_rounds", "2. Screening & workup burden", "Screening rounds offered (per person)", "person", 1,
  "screen_tp", "2. Screening & workup burden", "Test result: true positive (upper GI cancer, test +ve)", "count", 1,
  "screen_fp", "2. Screening & workup burden", "Test result: false positive (no upper GI cancer, test +ve)", "count", 0,
  "screen_fp_other_cancer", "2. Screening & workup burden", "  - of which another cancer present", "count", 1,
  "screen_tn", "2. Screening & workup burden", "Test result: true negative (no upper GI cancer, test -ve)", "count", 0,
  "screen_fn", "2. Screening & workup burden", "Test result: false negative (upper GI cancer, test -ve)", "count", 1,
  "fn_symptomatic_1y", "2. Screening & workup burden", "  - false negatives presenting with symptoms within a year", "count", 1,
  "fn_persistent", "2. Screening & workup burden", "  - missed cancers Galleri keeps missing (scenario p_fn_persist)", "count", 1,
  "investigated", "2. Screening & workup burden", "People investigated (EGD or more)", "count", 0,
  "egd_positive", "2. Screening & workup burden", "  - EGD +ve (cancer found at EGD)", "count", 1,
  "egd_negative_workup", "2. Screening & workup burden", "  - EGD -ve, sent on to further workup (Galleri) / EGD -ve (Kunzmann)", "count", 0,
  "egd", "2. Screening & workup burden", "EGDs", "count", 0,
  "ent", "2. Screening & workup burden", "ENT evaluations", "count", 0,
  "pet_ct", "2. Screening & workup burden", "PET CT scans", "count", 0,
  "eus", "2. Screening & workup burden", "EUS procedures", "count", 0,
  "eus_sampling", "2. Screening & workup burden", "EUS samplings", "count", 0,
  "tumor_board", "2. Screening & workup burden", "Tumor board discussions", "count", 0,
  "invasive_total", "2. Screening & workup burden", "Invasive procedures (EGD/EUS/sampling)", "count", 0,
  "mean_days", "2. Screening & workup burden", "Mean days to resolution if investigated", "ratio", 0,
  "cost_tests", "2. Screening & workup burden", "Screening, workup and complication cost per person ($, discounted)", "person", 2,
  "cost_tx", "2. Screening & workup burden", "Lifetime cancer care cost per person ($)", "person", 2,
  "cost_total", "2. Screening & workup burden", "Total lifetime cost per person ($)", "person", 2,
  "screen_detected", "3. Benefit", "Cancers found by screening/workup (all rounds)", "count", 1,
  "found_upper_gi", "3. Benefit", "  - upper GI (esophagus/stomach)", "count", 1,
  "found_other", "3. Benefit", "  - other sites", "count", 1,
  "clin_dx_lifetime", "3. Benefit", "Symptomatic diagnoses over lifetime (upper GI, net of other-site cancers found early)", "count", 1,
  "late_lifetime", "3. Benefit", "Late-stage diagnoses over lifetime", "count", 1,
  "late_averted", "3. Benefit", "Late-stage diagnoses averted vs reference", "count", 1,
  "ca_deaths", "3. Benefit", "Cancer deaths over lifetime", "count", 1,
  "ca_deaths_averted", "3. Benefit", "Cancer deaths averted vs reference", "count", 1,
  "oc_deaths", "3. Benefit", "Other-cause deaths (competing risk)", "count", 0,
  "ly_total", "3. Benefit", "Life years (undiscounted)", "count", 0,
  "ly_gained", "3. Benefit", "Life years gained vs reference", "count", 1,
  "overdiagnosed", "4. Harm", "Overdiagnosed cancers (screen-detected, would never have caused symptoms)", "count", 1,
  "overdiagnosed_ugi", "4. Harm", "  - upper GI", "count", 1,
  "overdiagnosed_other", "4. Harm", "  - other sites", "count", 1,
  "overdiagnosis_share", "4. Harm", "Share of screen-detected cancers overdiagnosed", "ratio", 3,
  "free_investigated", "4. Harm", "Cancer-free people investigated (false positives)", "count", 0,
  "egd_free", "4. Harm", "  - EGDs in cancer-free people", "count", 0,
  "ent_free", "4. Harm", "  - ENT evaluations in cancer-free people", "count", 0,
  "pet_ct_free", "4. Harm", "  - PET CTs in cancer-free people", "count", 0,
  "eus_free", "4. Harm", "  - EUS in cancer-free people", "count", 1,
  "eus_sampling_free", "4. Harm", "  - EUS samplings in cancer-free people", "count", 1,
  "incidental_followups_free", "4. Harm", "  - follow-ups of benign incidental PET findings", "count", 0,
  "perforations_free", "4. Harm", "  - EGD perforations in cancer-free people", "count", 2,
  "proc_deaths_free", "4. Harm", "  - deaths from complications in cancer-free people", "count", 3,
  "dose_sv_free", "4. Harm", "  - PET CT dose to cancer-free people (Sv)", "count", 2,
  "radiation_cancers_free", "4. Harm", "  - radiation-induced cancers in cancer-free people", "count", 3,
  "mean_days_fp", "4. Harm", "  - mean days in workup per cancer-free person", "ratio", 0,
  "days_free", "4. Harm", "  - person-days in workup, cancer-free people", "count", 0,
  "cost_fp", "4. Harm", "Workup and complication cost in cancer-free people ($ per person, discounted)", "person", 2,
  "qaly_fp", "4. Harm", "QALYs lost by cancer-free people worked up (procedures, complications, deaths, radiation)", "count", 2,
  "invasive_free", "4. Harm", "Invasive procedures in cancer-free people", "count", 0,
  "sae_total", "4. Harm", "Serious procedure complications (incl. perforation)", "count", 2,
  "perforations", "4. Harm", "  - EGD perforations", "count", 2,
  "perforation_deaths", "4. Harm", "  - fatal EGD perforations", "count", 3,
  "perforations_symptomatic", "4. Harm", "EGD perforations in the workup of symptomatic cancers", "count", 2,
  "incidental_followups", "4. Harm", "Follow-ups of incidental PET findings", "count", 0,
  "sae_free", "4. Harm", "  - in cancer-free people", "count", 2,
  "proc_deaths", "4. Harm", "Deaths from procedure complications", "count", 3,
  "dose_sv", "4. Harm", "Collective PET CT dose (Sv)", "count", 2,
  "radiation_cancers", "4. Harm", "Radiation-induced cancers (lifetime)", "count", 3,
  "missed_investigated", "4. Harm", "Cancers missed despite investigation (false reassurance)", "count", 1,
  "qaly_total", "5. QALYs & DALYs", "Lifetime QALYs (discounted)", "count", 0,
  "qaly_harms", "5. QALYs & DALYs", "QALYs lost to procedure harms", "count", 1,
  "qaly_net", "5. QALYs & DALYs", "Net QALYs gained vs reference", "count", 1,
  "daly_total", "5. QALYs & DALYs", "Cancer DALYs (YLL + YLD + harms)", "count", 0,
  "daly_harms", "5. QALYs & DALYs", "  - DALYs from procedure harms", "count", 1,
  "daly_averted", "5. QALYs & DALYs", "Net DALYs averted vs reference", "count", 1,
  "cost_per_daly", "6. Efficiency", "Cost per DALY averted vs reference ($)", "ratio", 0,
  "nn_per_ca_death_averted", "6. Efficiency", "People tested/worked up per cancer death averted", "ratio", 0,
  "investigated_per_found", "6. Efficiency", "People investigated per cancer found", "ratio", 1,
  "fp_per_found", "6. Efficiency", "False positives investigated per cancer found", "ratio", 1,
  "invasive_per_found", "6. Efficiency", "Invasive procedures per cancer found", "ratio", 1)

ALL_PARAMS <- c("sens_g_ugi", "sens_g_other", "fpr_g_es", "sens_egd", "p_egd_pos",
  "p_ent_route", "p_ent_pos", "p_pet_pb", "p_pet_incid", "p_eus_pos", "p_eus_neg", "p_miss",
  "p_incid_cancer", "p_late_screen", "or_late_clin", "sojourn_years", "m_early", "m_late", "p_perf_egd", "p_death_perf",
  "sae_ent", "sae_eus", "sae_eus_samp", "p_death_sae", "msv_pet", "risk_per_sv", "loss_rad_ca", "t_egd",
  "t_ent", "t_pet", "t_eus", "t_tb", "c_galleri", "c_kunz", "c_egd", "c_ent", "c_pet", "c_eus",
  "c_eus_samp", "c_tb", "c_symp_dx", "c_tx_first_early", "c_tx_first_late", "c_fu_early", "c_fu_late",
  "c_eol", "c_sae", "c_perf", "c_incid_fu", "sens_g_rep_ratio", "p_fn_persist", "galleri_cost_share", "u_mult_early", "u_mult_late", "d_proc", "d_sae", "dw_early", "dw_late")

keep_ordered <- function(q) {   # sampled values stay logically ordered
  q$u_mult_late <- min(q$u_mult_late, q$u_mult_early); q
}

KEY <- c("screen_detected", "late_averted", "ca_deaths_averted", "free_investigated", "overdiagnosed",
         "invasive_free", "sae_total", "qaly_net", "daly_averted")

make_analysis <- function(name, title, arms, owsa, cohort_n, unit, required, plots)
  ceakit::cea_analysis(
    name = name, title = title, arms = arms,
    run = function(p, data, settings) run_population(p, data, settings, arms)$total,
    cost = "cost_total", effect = "qaly_total", effect_label = "QALY",
    metrics = METRICS, derive = derive_galleri, required = required, cohort_n = cohort_n, unit = unit,
    owsa = owsa, key_metrics = KEY, constraints = keep_ordered, plots = plots)

UNIT <- "per 100,000 (base case: 50-year-old non-Hispanic White men)"
analyses <- list(
  galleri = make_analysis(
    "galleri", "Galleri: benefits and harms (annual screening)", c("Usual care", "Galleri"),
    c("Galleri", "Usual care"), 100000, UNIT, setdiff(ALL_PARAMS, c("c_kunz", "sens_egd")),
    list(workup_tree = plot_workup_tree, benefit_harm = plot_balance, icon_array = plot_icons,
         calibration = plot_calibration, markov_trace = plot_trace, cancer_states = plot_cancer_states,
         subgroups = plot_subgroups)),
  comparison = make_analysis(
    "comparison", "Usual care vs Kunzmann vs Galleri (annual screening)", ARMS,
    c("Galleri", "Kunzmann -> EGD"), 100000, UNIT, ALL_PARAMS,
    list(strategy_tree = plot_strategy_tree, benefit_harm = plot_balance, calibration = plot_calibration,
         markov_trace = plot_trace, cancer_states = plot_cancer_states, subgroups = plot_subgroups)))
