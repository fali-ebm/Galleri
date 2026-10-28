# Galleri model verification tests (run by run.R). Add one for every bug fixed.
proj <- ceakit::load_project(testthat::test_path("..", ".."))
E <- proj$env; p0 <- ceakit::param_list(proj$params); S <- proj$settings; D <- proj$data
ev <- function(name, p = p0) ceakit::evaluate(proj$analyses[[name]], p, D, S)

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

test_that("every cohort row's Markov trace conserves the cohort", {
  for (i in seq_len(nrow(D$cohort))) {
    curve <- E$baseline_for(D, D$cohort[i, ], S); mk <- E$markov_parts(curve, p0)
    expect_silent(ceakit::check_trace(ceakit::markov_trace(mk$n, mk$transition, c(1, 0, 0, 0, 0, 0) |>
      setNames(E$STATES))))
  }
})

test_that("regression: no further DALYs accrue after death (fixed double counting)", {
  curve <- E$baseline_for(D, D$cohort[1, ], S)
  V <- E$galleri_markov_values(curve, p0, S)
  expect_equal(unname(V[c("dead_cancer", "dead_other"), c("daly", "qaly", "cost", "ly")]), matrix(0, 2, 4))
})

test_that("no workup / no screening arms have no procedures or harms", {
  for (n in names(proj$analyses)) {
    m <- ev(n)
    expect_equal(unname(m[c("egd", "sae_total", "free_investigated", "screen_detected"), 1]), c(0, 0, 0, 0))
  }
})

test_that("a Galleri test with zero sensitivity and zero false positives only adds cost", {
  p <- p0; p$sens_g_ugi <- 0; p$sens_g_other <- 0; p$fpr_g_es <- 0
  m <- ev("screening", p)
  expect_equal(unname(m["qaly_total", "Galleri -> full workup"]), unname(m["qaly_total", 1]))
  expect_equal(unname(m["cost_total", "Galleri -> full workup"] - m["cost_total", 1]), p$c_galleri)
})

test_that("same total cancers in every arm (cancers are only moved between found and missed)", {
  for (n in names(proj$analyses)) {
    m <- ev(n)
    expect_equal(unname(m["screen_detected", ] + m["missed_total", ]), rep(unname(m["total_cancers", 1]), ncol(m)))
  }
})

test_that("worse late-stage survival increases cancer deaths", {
  p <- p0; p$m_late <- min(1, p0$m_late * 1.5)
  expect_gt(ev("galleri_positive", p)["ca_deaths", 1], ev("galleri_positive")["ca_deaths", 1])
})
