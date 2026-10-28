test_that("icer_frontier handles strong and extended dominance", {
  f <- icer_frontier(cost = c(0, 200, 300, 400), effect = c(0, 1, 2, 1.5),
                     strategies = c("A", "B", "C", "D"))
  st <- setNames(f$status, f$strategy)
  expect_equal(unname(st[c("A", "B", "C", "D")]),
               c("On frontier", "Extendedly dominated", "On frontier", "Dominated"))
  expect_equal(f$icer[f$strategy == "C"], 150)
})

test_that("EVPI is zero when one strategy always wins and positive otherwise", {
  set.seed(1)
  cost <- cbind(A = rep(0, 200), B = rep(100, 200))
  eff <- cbind(A = rep(0, 200), B = rep(1, 200))
  expect_equal(evpi(cost, eff, 1000)$evpi, 0)
  eff2 <- cbind(A = rep(0, 200), B = rnorm(200, 0.1, 1))
  expect_gt(evpi(cost, eff2, 1000)$evpi, 0)
  c2 <- ceac(cost, eff2, c(0, 1000))
  expect_equal(sum(c2$prob[c2$wtp == 1000]), 1)
})

test_that("sampling respects distributions and dirichlet groups", {
  P <- tidy_parameters(data.frame(
    parameter = c("p", "c", "d1", "d2", "f"), value = c(0.3, 100, 0.2, 0.5, 7),
    dist = c("beta", "gamma", "dirichlet", "dirichlet", "fixed"),
    se = c(NA, 10, NA, NA, NA), n = c(200, NA, 100, 100, NA), group = c(NA, NA, "g", "g", NA)))
  set.seed(2)
  s <- replicate(2000, unlist(sample_parameters(P)))
  expect_equal(mean(s["p", ]), 0.3, tolerance = 0.02)
  expect_equal(mean(s["c", ]), 100, tolerance = 0.02)
  expect_true(all(s["d1", ] + s["d2", ] <= 1))
  expect_true(all(s["f", ] == 7))
})

test_that("validation catches common input errors", {
  P <- tidy_parameters(data.frame(parameter = c("p_a", "c", "d1", "d2"), value = c(1.2, NA, 0.6, 0.6),
                                  dist = c("beta", "gamma", "dirichlet", "dirichlet"),
                                  n = c(10, NA, 10, 10), group = c(NA, NA, "g", "g"), status = "DUMMY"))
  v <- validate_parameters(P, required = c("p_a", "missing_one"))
  expect_true(any(grepl("missing rows", v$errors)))
  expect_true(any(grepl("outside 0-1", v$errors)))
  expect_true(any(grepl("no value", v$errors)))
  expect_true(any(grepl("sums to", v$errors)))
  expect_true(any(grepl("DUMMY", v$warnings)))
})

test_that("apply_counts converts counts to beta/dirichlet estimates", {
  P <- tidy_parameters(data.frame(parameter = c("p1", "q1", "q2"), value = 0.1, dist = "beta", n = 10))
  ct <- data.frame(branch = c("x", "x", "y", "y", "y"), outcome = letters[1:5],
                   parameter = c("p1", NA, "q1", "q2", NA), n = c(9, 1, 2, 3, 5))
  out <- apply_counts(P, ct)
  expect_equal(out$value[1], 9.5 / 11)
  expect_equal(out$dist[2:3], c("dirichlet", "dirichlet"))
  expect_error(apply_counts(P, transform(ct, n = c(9, NA, 2, 3, 5))), "partly filled")
})

test_that("run_cea runs end to end on a toy model", {
  P <- tidy_parameters(data.frame(parameter = c("c_tx", "eff"), value = c(1000, 0.05),
                                  dist = c("gamma", "normal"), se = c(200, 0.02)))
  a <- cea_analysis("toy", arms = c("Usual care", "New"),
    run = function(p, data, settings) {
      m <- rbind(cost = c(0, p$c_tx), qaly = c(1, 1 + p$eff)); colnames(m) <- c("Usual care", "New"); m },
    required = c("c_tx", "eff"))
  s <- default_settings(wtp = 50000, n_sim = 200)
  r <- run_cea(a, P, settings = s)
  expect_s3_class(r, "cea_result")
  expect_equal(unname(r$base["icer", "New"]), 1000 / 0.05)
  expect_equal(dim(r$psa)[3], 200)
  expect_true(all(c("c_tx", "eff") %in% r$evppi$parameter))
  expect_gt(r$evppi$evppi[r$evppi$parameter == "eff"], r$evppi$evppi[r$evppi$parameter == "c_tx"])
  sc <- run_scenarios(a, P, list(), s, data.frame(scenario = "cheap", parameter = "c_tx", value = 100))
  expect_equal(sc$icer[sc$scenario == "cheap" & sc$arm == "New"], 100 / 0.05)
  expect_equal(results_wide(r)$`Usual care`[results_wide(r)$Metric == "ICER vs Usual care ($ per QALY)"], "-")
})

# ---- Cross-validation against established packages (skipped if absent) ----
toy_result <- function(n_sim = 400) {
  P <- tidy_parameters(data.frame(parameter = c("c_tx", "eff", "noise"), value = c(1000, 0.02, 1),
                                  dist = c("gamma", "normal", "gamma"), se = c(400, 0.02, 0.2)))
  a <- cea_analysis("toy", arms = c("Usual care", "New A", "New B"),
    run = function(p, data, settings) {
      m <- rbind(cost = c(0, p$c_tx, 2 * p$c_tx * p$noise), qaly = c(1, 1 + p$eff, 1 + 1.6 * p$eff))
      colnames(m) <- c("Usual care", "New A", "New B"); m })
  run_cea(a, P, settings = default_settings(wtp = 50000, n_sim = n_sim, evppi = FALSE), owsa = FALSE)
}

test_that("ICERs, CEAC and EVPI match dampack and BCEA", {
  skip_if_not_installed("dampack"); skip_if_not_installed("BCEA")
  r <- toy_result(); m <- psa_matrices(r); w <- r$settings$wtp
  d <- dampack::calculate_icers(r$base["cost", ], r$base["qaly", ], r$arms)
  f <- icer_frontier(r$base["cost", ], r$base["qaly", ], r$arms)
  expect_equal(sort(d$ICER[!is.na(d$ICER)]), sort(f$icer[!is.na(f$icer)]))
  psa <- to_dampack(r)
  dc <- dampack::ceac(wtp = w, psa = psa)
  cc <- r$ceac[r$ceac$wtp == w, ]
  expect_equal(dc$Proportion[match(make.names(r$arms), dc$Strategy)], cc$prob[match(r$arms, cc$arm)])
  expect_equal(dampack::calc_evpi(psa, wtp = w)$EVPI, evpi(m$cost, m$effect, w)$evpi)
  b <- to_bcea(r, k = c(0, w))
  expect_equal(unname(b$evi[2]), evpi(m$cost, m$effect, w)$evpi)
})

test_that("built-in EVPPI agrees with voi", {
  skip_if_not_installed("voi")
  r <- toy_result(1000); m <- psa_matrices(r)
  ours <- evppi(m$cost, m$effect, m$parameters, 50000, method = "gam")
  theirs <- evppi(m$cost, m$effect, m$parameters, 50000, method = "voi")
  k <- match(ours$parameter, theirs$parameter)
  expect_equal(ours$evppi, theirs$evppi[k], tolerance = 0.02)
})
