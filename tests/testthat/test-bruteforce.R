test_that("bruteforce finds the obvious single breakpoint", {
  set.seed(6)
  n <- 400; bp_true <- 250; X <- matrix(1, n, 1)
  y <- c(rnorm(bp_true, 0, 1), rnorm(n - bp_true, 8, 1))
  res <- bruteforce(y, X, gaussian_homo())
  expect_true(abs(res$bp - bp_true) <= 3)
})
