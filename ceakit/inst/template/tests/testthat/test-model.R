# Project verification tests (run by run.R). Add a test for every bug you fix.
proj <- ceakit::load_project(testthat::test_path("..", ".."))
p0 <- ceakit::param_list(proj$params)
a <- proj$analyses$base
ev <- function(p) ceakit::evaluate(a, p, proj$data, proj$settings)

test_that("inputs are valid", {
  expect_length(ceakit::check_project(proj)$errors, 0)
})

test_that("a useless test (sensitivity 0, specificity 1) only adds cost", {
  p <- p0; p$sens <- 0; p$spec <- 1
  m <- ev(p)
  expect_equal(unname(m["qaly", "Screening"]), unname(m["qaly", "Usual care"]))
  expect_gt(m["cost", "Screening"], m["cost", "Usual care"])
})

test_that("more prevalence means more cases detected", {
  p <- p0; p$p_prev <- p0$p_prev * 2
  expect_gt(ev(p)["detected", "Screening"], ev(p0)["detected", "Screening"])
})

test_that("Markov trace conserves the cohort", {
  mk <- proj$env$markov_parts(p0, proj$data, proj$settings)
  tr <- ceakit::markov_trace(mk$n, mk$transition, c(well = 0.5, undiagnosed = 0.5, treated = 0, dead = 0))
  expect_silent(ceakit::check_trace(tr))
})
