test_that("glmseg one-breakpoint result matches brute force", {
  set.seed(9)
  n <- 500; bp_true <- 300; X <- matrix(1, n, 1)
  y <- c(rnorm(bp_true, 1, 1), rnorm(n - bp_true, 7, 1))
  fam <- gaussian_homo()
  fit <- glmseg(y, X, K = 2, family = fam, init = "binseg")
  bf <- bruteforce(y, X, fam)
  expect_true(abs(fit$bp - bf$bp) <= 3)
  expect_equal(fit$loglik, bf$loglik, tolerance = 1e-4)
})

test_that("glmseg builds default intercept X when X is NULL", {
  set.seed(10)
  y <- c(rnorm(150, 0, 1), rnorm(150, 5, 1))
  fit <- glmseg(y, K = 2, family = gaussian_homo(), init = "binseg")
  expect_equal(fit$p, 1L)
  expect_length(fit$bp, 1L)
})

test_that("glmseg handles the null model K = 1 (no breakpoints)", {
  set.seed(21)
  y <- rnorm(300, 2, 1); X <- matrix(1, 300, 1)
  fit <- glmseg(y, X, K = 1, family = gaussian_homo(), init = "binseg")
  expect_length(fit$bp, 0L)
  expect_equal(fit$K, 1L)
  expect_equal(fit$loglik,
               segmentation_loglik(gaussian_homo(), y, X, integer(0)),
               tolerance = 1e-8)
})

test_that("glmseg works end-to-end with the default fused_lasso initialization", {
  set.seed(31)
  n <- 500; bp_true <- 300; X <- matrix(1, n, 1)
  y <- c(rnorm(bp_true, 0, 1), rnorm(n - bp_true, 9, 1))
  fit <- glmseg(y, X, K = 2, family = gaussian_homo(), init = "fused_lasso")
  expect_length(fit$bp, 1L)
  expect_true(abs(fit$bp - bp_true) <= 10)
})
