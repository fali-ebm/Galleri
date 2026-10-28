test_that("tree_paths multiplies probabilities and sums node attributes", {
  e <- data.frame(from = c("root", "root", "a", "a"), to = c("a", "b", "a1", "a2"),
                  prob = c(0.3, 0.7, 0.4, 0.6))
  n <- data.frame(id = c("a", "a1", "b"), cost = c(10, 5, 1))
  d <- tree_paths(e, n)
  expect_equal(sum(d$prob), 1)
  expect_equal(d$prob[d$leaf == "a1"], 0.12)
  expect_equal(d$cost[d$leaf == "a1"], 15)
  expect_equal(d$cost[d$leaf == "a2"], 10)
  expect_error(tree_paths(transform(e, prob = c(0.3, 0.6, 0.4, 0.6))), "do not sum to 1")
})

test_that("markov_values matches the closed-form answer", {
  P <- matrix(c(0.9, 0.1, 0, 1), 2, byrow = TRUE, dimnames = list(c("alive", "dead"), c("alive", "dead")))
  R <- matrix(c(1, 0), 2, 1, dimnames = list(c("alive", "dead"), "ly"))
  v <- markov_values(10, function(t) P, function(t, P) R)
  expect_equal(v["alive", "ly"], sum(0.9^(0:9)))
  vd <- markov_values(10, function(t) P, function(t, P) R, disc = c(ly = 0.05))
  expect_equal(vd["alive", "ly"], sum((0.9 / 1.05)^(0:9)))
  tr <- markov_trace(10, function(t) P, c(alive = 1, dead = 0))
  expect_silent(check_trace(tr))
  expect_equal(tr$alive[11], 0.9^10)
  bad <- P; bad[1, 1] <- 0.95
  expect_error(markov_values(2, function(t) bad, function(t, P) R), "do not sum to 1")
})

test_that("alive_fraction applies half-cycle correction", {
  P <- matrix(c(0.8, 0.2, 0, 1), 2, byrow = TRUE, dimnames = list(c("a", "d"), c("a", "d")))
  expect_equal(alive_fraction(P, "d"), c(0.9, 0))
})

test_that("lookup_by_age falls back to 'all' and checks coverage", {
  tab <- data.frame(group = c("all", "all", "hi"), sex = "all", age_from = c(0, 50, 0),
                    age_to = c(49, 100, 100), q = c(0.01, 0.02, 0.5))
  expect_equal(lookup_by_age(tab, list(group = "lo", sex = "male"), c(40, 60))$q, c(0.01, 0.02))
  expect_equal(lookup_by_age(tab, list(group = "hi", sex = "male"), 60)$q, 0.5)
  expect_error(lookup_by_age(tab, list(group = "hi", sex = "male"), 101), "do not cover")
})
