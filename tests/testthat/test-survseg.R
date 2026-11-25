test_that("survseg recovers a breakpoint (binseg and fused_lasso)", {
  skip_if_not_installed("survival")
  set.seed(6)
  n <- 600; bp_true <- 300; X <- matrix(1, n, 1)
  d <- sim_survseg(X, par = c(0.3, 3.0, 1.3, 1.3), bp = bp_true,
                   cens_scale = c(80, 80))
  fit_bs <- survseg(d$time, d$status, X, K = 2, init = "binseg")
  expect_s3_class(fit_bs, "survseg")
  expect_length(fit_bs$bp, 1L)
  expect_true(abs(fit_bs$bp - bp_true) <= 60)
  fit_fl <- survseg(d$time, d$status, X, K = 2, init = "fused_lasso")
  expect_length(fit_fl$bp, 1L)
  expect_true(abs(fit_fl$bp - bp_true) <= 80)
})

test_that("survseg handles the null model K=1 and select_k recovers K", {
  skip_if_not_installed("survival")
  set.seed(7)
  n <- 900; X <- matrix(1, n, 1)
  d <- sim_survseg(X, par = c(0.2, 2.6, 0.2, 1.3, 1.3, 1.3),
                   bp = c(300, 600), cens_scale = c(90, 90, 90))
  f1 <- survseg(d$time, d$status, X, K = 1, init = "binseg")
  expect_length(f1$bp, 0L)
  res <- survseg_select_k(d$time, d$status, X, K_max = 3, init = "binseg")
  expect_equal(res$K, 3L)
})

test_that("survseg_bic mbic adds the (K-1)*log(n) breakpoint penalty", {
  skip_if_not_installed("survival")
  set.seed(6)
  n <- 500; bp_true <- 250; X <- matrix(1, n, 1)
  d <- sim_survseg(X, par = c(0.5, 3.0, 1.3, 1.3), bp = bp_true,
                   cens_scale = c(80, 80))
  fit <- survseg(d$time, d$status, X, K = 2, init = "binseg")
  expect_equal(survseg_bic(fit, "mbic"),
               survseg_bic(fit, "bic") + (fit$K - 1L) * log(fit$n))
  expect_equal(survseg_bic(fit), survseg_bic(fit, "bic"))
})

test_that("survseg_select_k with mbic recovers the true K", {
  skip_if_not_installed("survival")
  set.seed(7)
  n <- 900; X <- matrix(1, n, 1)
  d <- sim_survseg(X, par = c(0.2, 2.6, 0.2, 1.3, 1.3, 1.3),
                   bp = c(300, 600), cens_scale = c(90, 90, 90))
  res <- survseg_select_k(d$time, d$status, X, K_max = 3, init = "binseg",
                          penalty = "mbic")
  expect_equal(res$K, 3L)
})
