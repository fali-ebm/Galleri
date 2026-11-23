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
  C <- E$subgroup_cohort(60, "male", "nhw", S); curve <- E$baseline_for(D, C[1, ], S)
  t_hi <- E$target_ugi(curve, C[1, ], S); t_lo <- E$target_ugi(curve, C[2, ], S)
  expect_equal(S$kunz_share * t_hi + (1 - S$kunz_share) * t_lo, curve$inc_eso + curve$inc_stom)
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
  p <- p0; p$sens_g_ugi <- 1                                               # perfect sensitivity: no false negatives
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
