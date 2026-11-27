# Galleri model verification tests (run by run.R). Add one for every bug fixed.
proj <- ceakit::load_project(testthat::test_path("..", ".."))
E <- proj$env; p0 <- ceakit::param_list(proj$params); S <- proj$settings; D <- proj$data
ev <- function(name, p = p0, s = S) ceakit::evaluate(proj$analyses[[name]], p, D, s)
row1 <- D$cohort[1, ]
sim <- function(arm, row = row1, p = p0, s = S) {
  prep <- E$prepare_row(row, D, p, s)
  E$simulate_arm(arm, row, prep$curve, prep$cal, p, s, E$workup_summary(p), prep$V)
}

test_that("inputs are valid", {
  expect_length(ceakit::check_project(proj)$errors, 0)
})

test_that("workup tree probabilities sum to 1 and procedures are counted per path", {
  d <- E$workup_paths(p0)
  expect_equal(sum(d$prob), 1)
  expect_equal(d$n_egd, rep(1, nrow(d)))
  expect_equal(d$n_samp[d$leaf == "eus_pos"], 1)
  expect_equal(d$n_tb[d$leaf == "egd_pos"], 0)
})

test_that("calibration: usual care reproduces SEER incidence at every age, for every subgroup", {
  for (race in S$subgroups$races) for (sex in S$subgroups$sexes) for (age in c(50, 70, 90)) {
    C <- E$subgroup_cohort(age, sex, race, S)
    ct <- E$calibration_table(p0, D, S, "Usual care", C)
    expect_equal(ct$model_clinical, ct$seer, tolerance = 1e-8, label = paste(race, sex, age))
    expect_equal(unique(ct$clipped), 0)
  }
})

test_that("stage calibration: symptomatic late-stage share matches SEER by site mix, race and sex", {
  for (race in c("nhw", "nhb", "hisp")) for (sex in S$subgroups$sexes) {
    C <- E$subgroup_cohort(55, sex, race, S)
    ct <- E$calibration_table(p0, D, S, "Usual care", C)
    expect_equal(ct$model_late_share, ct$seer_late_share, tolerance = 1e-10)
    curve <- E$baseline_for(D, C[1, ], S)
    expect_true(all(ct$seer_late_share >= pmin(curve$late_eso, curve$late_stom) - 1e-12 &
                    ct$seer_late_share <= pmax(curve$late_eso, curve$late_stom) + 1e-12))
  }
  p <- p0; p$or_late_clin <- 2
  expect_gt(mean(E$calibration_table(p, D, S, "Usual care")$model_late_share),
            mean(E$calibration_table(p0, D, S, "Usual care")$model_late_share))
})

test_that("calibration holds when natural-history inputs change (PSA draws re-calibrate)", {
  p <- p0; p$sojourn_years <- 4; p$m_late <- 0.5
  ct <- E$calibration_table(p, D, S, "Usual care")
  expect_equal(ct$model_clinical, ct$seer, tolerance = 1e-8)
})

test_that("Kunzmann groups average back to the subgroup's SEER incidence", {
  for (sx in c("male", "female")) {
    C <- E$subgroup_cohort(60, sx, "nhw", S); curve <- E$baseline_for(D, C[1, ], S)
    t_hi <- E$target_ugi(curve, C[1, ], S); t_lo <- E$target_ugi(curve, C[2, ], S)
    k <- S[[paste0("kunz_share_", sx)]]
    expect_equal(C$share, c(k, 1 - k))
    expect_equal(k * t_hi + (1 - k) * t_lo, curve$inc_eso + curve$inc_stom)
  }
  expect_gte(S$kunz_share_male, S$kunz_share_female)
})

test_that("every arm's trace conserves the cohort", {
  for (a in E$ARMS) {
    tr <- sim(a)$trace
    expect_equal(rowSums(tr), rep(1, nrow(tr)), tolerance = 1e-10)
    expect_true(all(tr >= -1e-12))
  }
})

test_that("screening stops at the stop age and repeats at the interval", {
  s <- S; s$screen_stop_age <- 60
  expect_equal(unname(sim("Galleri", s = s)$ev[["rounds"]]), 11)           # ages 50-60
  s$screen_interval <- 2
  expect_equal(unname(sim("Galleri", s = s)$ev[["rounds"]]), 6)            # 50, 52, ..., 60
  row90 <- E$subgroup_cohort(90, "male", "nhw", S)[1, ]
  expect_equal(unname(sim("Galleri", row = row90)$ev[["rounds"]]), 0)      # no screening after 85
})

test_that("regression: no lead-time penalty - detection with no stage shift changes no deaths", {
  D2 <- D; D2$baseline_risks$late_eso <- 0.6; D2$baseline_risks$late_stom <- 0.6   # constant clinical mix
  p <- p0; p$p_late_screen <- 0.6; p$sens_g_other <- 0
  m <- ceakit::evaluate(proj$analyses$galleri, p, D2, S)
  expect_equal(unname(m["ca_deaths", "Galleri"]), unname(m["ca_deaths", "Usual care"]), tolerance = 1e-6)
})

test_that("stage shift at screening reduces cancer deaths", {
  m <- ev("galleri")
  expect_lt(m["ca_deaths", "Galleri"], m["ca_deaths", "Usual care"])
})

test_that("regression: no further DALYs or QALYs accrue after death", {
  prep <- E$prepare_row(row1, D, p0, S)
  for (V in prep$V) expect_equal(unname(V[c("dead_cancer", "dead_other"), c("daly", "qaly", "cost", "ly")]), matrix(0, 2, 4))
})

test_that("the usual care arm has no procedures or harms", {
  for (n in names(proj$analyses)) {
    m <- ev(n)
    expect_equal(unname(m[c("egd", "sae_total", "free_investigated", "screen_detected", "galleri_tests"), 1]),
                 c(0, 0, 0, 0, 0))
  }
})

test_that("a Galleri test with zero sensitivity and zero false positives only adds test costs", {
  p <- p0; p$sens_g_ugi <- 0; p$sens_g_other <- 0; p$fpr_g_es <- 0
  m <- ev("galleri", p)
  expect_equal(unname(m["qaly_total", "Galleri"]), unname(m["qaly_total", 1]))
  expect_equal(unname(m["ca_deaths", "Galleri"]), unname(m["ca_deaths", 1]))
  expect_gt(m["cost_total", "Galleri"], m["cost_total", 1])
  # one-off test: the extra cost is exactly the test price
  s <- S; s$screen_stop_age <- 50
  m1 <- ev("galleri", p, s)
  expect_equal(unname(m1["cost_total", "Galleri"] - m1["cost_total", 1]), p$c_galleri)
})

test_that("more screening rounds find more cancers", {
  s <- S; s$screen_stop_age <- 50
  expect_gt(ev("galleri")["screen_detected", "Galleri"], ev("galleri", s = s)["screen_detected", "Galleri"])
})

test_that("worse late-stage survival increases cancer deaths", {
  p <- p0; p$m_late <- min(1, p0$m_late * 1.5)
  expect_gt(ev("galleri", p)["ca_deaths", 1], ev("galleri")["ca_deaths", 1])
})

test_that("both analyses give the same usual care and Galleri results", {
  a <- ev("galleri"); b <- ev("comparison")
  expect_equal(a[, c("Usual care", "Galleri")], b[, c("Usual care", "Galleri")])
})

test_that("the Galleri arm runs each positive through the workup flowchart", {
  m <- ev("galleri")[, "Galleri"]; d <- E$workup_paths(p0)
  expect_equal(m[["egd"]], m[["investigated"]])                       # every positive has an EGD
  s_egd <- 1 - p0$p_perf_egd * p0$p_death_perf                          # a fatal perforation ends the workup
  expect_equal(m[["pet_ct"]], m[["investigated"]] * s_egd * sum(d$prob * d$n_pet))
  p <- p0; p$p_miss <- 0.5                                              # follow-up counts matter
  expect_lt(ev("galleri", p)["screen_detected", "Galleri"], m[["screen_detected"]])
  expect_equal(m[["screen_detected"]] + m[["missed_investigated"]] + m[["free_investigated"]],
               m[["investigated"]], tolerance = 1e-12)
})

test_that("subgroup grid covers every age x sex x race/ethnicity and screening stops at 85", {
  s <- S; s$subgroups <- list(ages = c(50, 90), sexes = c("male", "female"), races = c("nhw", "hisp"))
  g <- E$subgroup_grid(p0, D, s, E$ARMS)
  expect_equal(nrow(g), 2 * 2 * 2 * length(E$ARMS))
  expect_true(all(g$cancers_found[g$start_age == 90] == 0))
})

test_that("Galleri 2x2: every test is exactly one of TP, FP, TN, FN, and TP + FP are the people investigated", {
  for (arm in c("Galleri", "Kunzmann -> EGD")) {
    m <- ev("comparison")[, arm]
    tested <- if (arm == "Galleri") m[["galleri_tests"]] else m[["investigated"]]
    expect_equal(m[["screen_tp"]] + m[["screen_fp"]] + m[["screen_tn"]] + m[["screen_fn"]], tested, tolerance = 1e-10)
    if (arm == "Galleri") expect_equal(m[["screen_tp"]] + m[["screen_fp"]], m[["investigated"]], tolerance = 1e-10)
    expect_lte(m[["fn_symptomatic_1y"]], m[["screen_fn"]])
  }
  p <- p0; p$sens_g_ugi <- 1; p$sens_g_rep_ratio <- 1                     # perfect sensitivity: no false negatives
  expect_equal(unname(ev("galleri", p)["screen_fn", "Galleri"]), 0)
  p <- p0; p$sens_g_other <- 0; p$fpr_g_es <- 0                            # perfect specificity: no false positives
  expect_equal(unname(ev("galleri", p)["screen_fp", "Galleri"]), 0)
  m <- ev("galleri")[, "Galleri"]                                          # EGD +ve / -ve split of those investigated
  expect_equal(m[["egd_positive"]] + m[["egd_negative_workup"]],
               m[["investigated"]] * (1 - p0$p_perf_egd * p0$p_death_perf), tolerance = 1e-10)
})

test_that("EGD perforation: counted on every EGD, fatal share follows p_death_perf, and it costs QALYs and money", {
  for (arm in c("Galleri", "Kunzmann -> EGD")) {
    m <- ev("comparison")[, arm]
    expect_equal(m[["perforations"]], m[["egd"]] * p0$p_perf_egd, tolerance = 1e-10)
    expect_equal(m[["perforation_deaths"]], m[["perforations"]] * p0$p_death_perf, tolerance = 1e-10)
  }
  p <- p0; p$p_perf_egd <- 0
  m0 <- ev("comparison", p); m1 <- ev("comparison")
  expect_equal(unname(m0["perforations", "Kunzmann -> EGD"]), 0)
  expect_gt(m0["qaly_total", "Kunzmann -> EGD"], m1["qaly_total", "Kunzmann -> EGD"])
  expect_lt(m0["cost_total", "Kunzmann -> EGD"], m1["cost_total", "Kunzmann -> EGD"])
  p <- p0; p$p_death_perf <- 0.5                                           # deadlier perforation -> more procedure deaths
  expect_gt(ev("comparison", p)["proc_deaths", "Kunzmann -> EGD"], m1["proc_deaths", "Kunzmann -> EGD"])
})

test_that("overdiagnosis equals the excess lifetime diagnoses vs usual care (no procedure deaths)", {
  # clin_dx_lifetime in a screening arm is net of the symptomatic diagnoses that
  # other-site cancers found by the workup would have had, so all finds count
  p <- p0; p$p_perf_egd <- 0; p$p_death_sae <- 0           # procedure deaths would change later onset
  m <- ev("comparison", p)
  for (arm in c("Galleri", "Kunzmann -> EGD")) {
    excess <- m["found_upper_gi", arm] + m["found_other", arm] + m["clin_dx_lifetime", arm] -
              m["clin_dx_lifetime", "Usual care"]
    expect_equal(unname(m["overdiagnosed", arm]), unname(excess), tolerance = 1e-9)
  }
  m <- ev("galleri")
  expect_gt(m["overdiagnosed", "Galleri"], 0)
  expect_lt(m["overdiagnosed", "Galleri"], m["screen_detected", "Galleri"])
  expect_equal(unname(m["overdiagnosed", 1]), 0)
  # older people and slower-growing cancers are overdiagnosed more often
  older <- row1; older$start_age <- 75
  share <- function(r, p = p0) { s <- sim("Galleri", r, p); (s$ev[["overdx_ugi"]] + s$ev[["overdx_other"]]) /
                                   (s$ev[["found_ugi"]] + s$ev[["found_other"]]) }
  expect_gt(share(older), share(row1))
  p <- p0; p$sojourn_years <- p0$sojourn_years * 2
  expect_gt(share(row1, p), share(row1))
})

test_that("false-positive burden: counted only in cancer-free people and consistent with the totals", {
  m <- ev("galleri")[, "Galleri"]
  expect_equal(m[["egd_free"]], m[["free_investigated"]])            # every false positive has an EGD
  for (k in c("ent", "pet_ct", "eus", "eus_sampling", "perforations"))
    expect_lte(m[[paste0(k, "_free")]], m[[k]] + 1e-15)
  expect_lte(m[["proc_deaths_free"]], m[["proc_deaths"]] + 1e-15)
  expect_gt(m[["qaly_fp"]], 0); expect_gt(m[["cost_fp"]], 0)
  p <- p0; p$fpr_g_es <- 0                                            # no false positives -> no burden
  m0 <- ev("galleri", p)[, "Galleri"]
  for (k in c("free_investigated", "egd_free", "pet_ct_free", "perforations_free", "qaly_fp", "cost_fp"))
    expect_equal(m0[[k]], 0)
  p <- p0; p$c_perf <- p0$c_perf * 2                                   # perforation cost feeds the FP cost
  expect_gt(ev("galleri", p)["cost_fp", "Galleri"], m[["cost_fp"]])
})

test_that("cause-deleted other-cause mortality removes SEER upper GI cancer deaths from the life table", {
  B <- D$baseline_risks
  expect_true(all(B$mort_other < B$mort_all))
  expect_equal(B$mort_other, 1 - (1 - B$mort_all) / (1 - B$mort_ugi), tolerance = 1e-9)
})

test_that("costs are inflated to the cost year with CPI medical care", {
  P <- read.csv(testthat::test_path("..", "..", "inputs", "parameters.csv"))
  cpi <- read.csv(testthat::test_path("..", "..", "data-raw", "raw", "cpi_medical_care.csv"))
  r <- P[!is.na(P$price_year), ]
  f <- cpi$index[match(S$cost_year, cpi$year)] / cpi$index[match(pmin(r$price_year, S$cost_year), cpi$year)]
  expect_equal(r$value, round(r$price_value * f))
  expect_gt(nrow(r), 10)
})

test_that("repeat-round Galleri sensitivity: lower ratio finds fewer cancers; no effect with one round", {
  p <- p0; p$sens_g_rep_ratio <- 1
  expect_gt(ev("galleri", p)["screen_tp", "Galleri"], ev("galleri")["screen_tp", "Galleri"])
  s1 <- S; s1$screen_stop_age <- 50
  expect_equal(ev("galleri", p, s1)["screen_tp", "Galleri"], ev("galleri", s = s1)["screen_tp", "Galleri"])
})

test_that("scenario: missed cancers that stay missed are never found by later Galleri rounds", {
  p <- p0; p$p_fn_persist <- 1
  s <- sim("Galleri", p = p); s0 <- sim("Galleri")
  expect_gt(s$ev[["fn_persist"]], 0); expect_equal(s0$ev[["fn_persist"]], 0)
  expect_equal(unname(rowSums(s$trace)), rep(1, nrow(s$trace)), tolerance = 1e-12)
  expect_lt(s$ev[["found_ugi"]], s0$ev[["found_ugi"]])
  expect_equal(unname(s0$trace[, "pre_m"]), rep(0, nrow(s0$trace)))
  m <- ev("galleri", p)[, "Galleri"]                                       # 2x2 still adds up
  expect_equal(m[["screen_tp"]] + m[["screen_fp"]] + m[["screen_tn"]] + m[["screen_fn"]], m[["galleri_tests"]], tolerance = 1e-10)
  # usual care is unchanged, so calibration still holds
  expect_equal(ev("galleri", p)[, "Usual care"], ev("galleri")[, "Usual care"])
})

test_that("scenario: Galleri price sharing scales only the test cost", {
  s1 <- S; s1$screen_stop_age <- 50                                         # one round: everyone tested once
  p <- p0; p$galleri_cost_share <- 0.5
  d <- ev("galleri", s = s1)["cost_total", "Galleri"] - ev("galleri", p, s1)["cost_total", "Galleri"]
  expect_equal(unname(d), 0.5 * p0$c_galleri, tolerance = 1e-9)
})

test_that("incidental PET findings are followed up, and benign ones add false-positive burden", {
  m <- ev("galleri")[, "Galleri"]; d <- E$workup_paths(p0); s_egd <- 1 - p0$p_perf_egd * p0$p_death_perf
  expect_equal(m[["incidental_followups"]], m[["investigated"]] * s_egd * sum(d$prob * d$n_incid), tolerance = 1e-10)
  expect_gt(m[["incidental_followups_free"]], 0)
  p <- p0; p$c_incid_fu <- p0$c_incid_fu * 2
  expect_gt(ev("galleri", p)["cost_fp", "Galleri"], m[["cost_fp"]])
})

test_that("symptomatic diagnoses carry EGD perforation risk in every arm", {
  m <- ev("comparison")
  expect_equal(unname(m["perforations_symptomatic", "Usual care"]),
               unname(m["clin_dx_lifetime", "Usual care"]) * p0$p_perf_egd, tolerance = 1e-10)
  p <- p0; p$p_perf_egd <- 0
  expect_gt(ev("comparison", p)["qaly_total", "Usual care"], m["qaly_total", "Usual care"])
})

test_that("workbook importer: the committed workbook validates, and bad edits are rejected without writing", {
  xl <- testthat::test_path("..", "..", "Galleri_model_inputs.xlsx")
  skip_if_not(file.exists(xl), "no inputs workbook (python3 tools/export_inputs_xlsx.py)")
  skip_if_not_installed("openxlsx")
  source(testthat::test_path("..", "..", "tools", "import_inputs.R"), local = TRUE)
  proj_dir <- testthat::test_path("..", "..")
  r <- import_workbook(xl, proj_dir, apply = FALSE, quiet = TRUE)
  expect_true(r$ok); expect_length(r$errors, 0)
  if (length(r$changes)) skip(paste("workbook is out of date with the CSVs - re-export it:", r$changes[1]))
  wb <- openxlsx::loadWorkbook(xl)
  P <- openxlsx::read.xlsx(wb, "2 EGD"); i <- which(P$Parameter == "p_egd_pos")
  openxlsx::writeData(wb, "2 EGD", 1.5, startCol = which(names(P) == "Value"), startRow = i + 1)
  bad <- tempfile(fileext = ".xlsx"); openxlsx::saveWorkbook(wb, bad)
  before <- tools::md5sum(file.path(proj_dir, "inputs", "parameters.csv"))
  r2 <- import_workbook(bad, proj_dir, apply = TRUE, quiet = TRUE)
  expect_false(r2$ok)
  expect_true(any(grepl("p_egd_pos: proportions must be 0-1", r2$errors)))
  expect_equal(tools::md5sum(file.path(proj_dir, "inputs", "parameters.csv")), before)
})

test_that("false-positive burden: components add up, pathway TTO replaces per-procedure loss, scenario 0 restores it", {
  m <- ev("galleri")[, "Galleri"]
  expect_equal(m[["qaly_fp"]], m[["fp_qaly_pathway"]] + m[["fp_qaly_worry"]] + m[["fp_qaly_complications"]] +
                 m[["fp_qaly_deaths"]] + m[["fp_qaly_radiation"]], tolerance = 1e-10)
  expect_equal(m[["daly_fp"]], m[["fp_daly_anxiety"]] + m[["fp_daly_worry"]] + m[["fp_daly_complications"]] +
                 m[["fp_daly_deaths"]] + m[["fp_daly_radiation"]], tolerance = 1e-10)
  expect_equal(m[["fp_worry_cases"]], m[["free_investigated"]] * p0$p_worry, tolerance = 1e-10)
  # pathway loss per false positive lies between the smallest and largest Matza values
  per <- m[["fp_qaly_pathway"]] / m[["free_investigated"]]
  expect_gt(per, 0.6 * p0$dfp_endo); expect_lt(per, p0$dfp_ext)        # discounted over the screening years
  # switching to per-procedure loss: pathway QALYs = procedures x d_proc (discounted), harms fall
  p <- p0; p$fp_tto <- 0; m0 <- ev("galleri", p)[, "Galleri"]
  n_proc <- m0[["egd_free"]] + m0[["ent_free"]] + m0[["pet_ct_free"]] + m0[["eus_free"]] + m0[["incidental_followups_free"]]
  expect_lt(m0[["fp_qaly_pathway"]], n_proc * p0$d_proc)        # discounted
  expect_gt(m0[["fp_qaly_pathway"]], 0.4 * n_proc * p0$d_proc)
  expect_lt(m0[["qaly_harms"]], m[["qaly_harms"]])
  # more worry or a heavier anxiety weight only add burden
  p <- p0; p$dw_anx <- 0.133; expect_gt(ev("galleri", p)["daly_fp", "Galleri"], m[["daly_fp"]])
  p <- p0; p$p_worry <- 0; expect_equal(unname(ev("galleri", p)["fp_qaly_worry", "Galleri"]), 0)
  # no false positives, no false-positive burden
  p <- p0; p$fpr_g_es <- 0; expect_equal(unname(ev("galleri", p)["daly_fp", "Galleri"]), 0)
})

test_that("accuracy by age, sex and race: 2x2 adds up and PPV follows SEER risk", {
  a <- E$accuracy_table(p0, D, S)
  expect_equal(a$tp + a$fp + a$tn + a$fn, a$tests, tolerance = 1e-8)
  expect_equal(nrow(a), length(S$subgroups$sexes) * length(S$subgroups$races) * (S$screen_stop_age - 50 + 1))
  # sensitivity is a test property: round 1 vs later rounds, the same in every group
  r1 <- a$sensitivity[a$age == 50]; r2 <- a$sensitivity[a$age == 70]
  expect_equal(r1, rep(p0$sens_g_ugi, length(r1)), tolerance = 1e-8)
  expect_equal(r2, rep(p0$sens_g_ugi * p0$sens_g_rep_ratio, length(r2)), tolerance = 1e-8)
  # PPV is higher where risk is higher: men vs women, age 80 vs 55
  g <- function(sx, ag) a$ppv_upper_gi[a$sex == sx & a$race_eth == "all" & a$age == ag]
  expect_gt(g("male", 70), g("female", 70)); expect_gt(g("male", 80), g("male", 55))
})

test_that("Kunzmann arm: one-time score and one-time EGD for the flagged share", {
  m <- ev("comparison")
  expect_equal(unname(m["kunz_scores", "Kunzmann -> EGD"]), 1, tolerance = 1e-9)
  expect_equal(unname(m["screen_rounds", "Kunzmann -> EGD"]), S$kunz_share_male, tolerance = 1e-6)
  S2 <- S; S2$kunz_egd_rounds <- 3
  expect_equal(unname(ev("comparison", s = S2)["screen_rounds", "Kunzmann -> EGD"]), 3 * S$kunz_share_male, tolerance = 1e-3)
})
