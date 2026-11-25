test_that("fused-lasso proposes candidates near a strong breakpoint", {
  set.seed(7)
  n <- 500; bp_true <- 300
  y <- c(rnorm(bp_true, 0, 1), rnorm(n - bp_true, 10, 1))
  cand <- init_fused_lasso(y, K = 2)
  expect_true(length(cand) >= 1)
  # Fused-Lasso is a coarse initialization later refined by max-EM; thinning
  # (almost_unique_tol) keeps the left-most candidate of a cluster, so a match
  # within a few dozen points of the true breakpoint is the meaningful contract.
  expect_true(min(abs(cand - bp_true)) <= 50)
})

test_that("init_fused_lasso accepts an explicit glmnet_family and still recovers gaussian breakpoints", {
  set.seed(7)
  n <- 500; bp_true <- 300
  y <- c(rnorm(bp_true, 0, 1), rnorm(n - bp_true, 10, 1))
  cand <- init_fused_lasso(y, K = 2, glmnet_family = "gaussian")
  expect_true(length(cand) >= 1)
  expect_true(min(abs(cand - bp_true)) <= 50)
})
