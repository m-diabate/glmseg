test_that("survseg S3 accessors and aft_to_ph are coherent", {
  skip_if_not_installed("survival")
  set.seed(8)
  n <- 500; bp_true <- 250; X <- matrix(1, n, 1)
  d <- sim_survseg(X, par = c(0.5, 3.0, 1.3, 1.3), bp = bp_true,
                   cens_scale = c(80, 80))
  fit <- survseg(d$time, d$status, X, K = 2, init = "binseg")
  expect_equal(as.numeric(logLik(fit)), fit$loglik)
  expect_equal(attr(logLik(fit), "df"), (fit$p + 1L) * fit$K)
  expect_equal(BIC(fit), -2 * fit$loglik + (fit$p + 1L) * fit$K * log(n))
  ph <- aft_to_ph(fit)
  # PH coefficient = -shape * AFT coefficient
  expect_equal(ph[[1]], -fit$shape[1] * fit$betas[[1]], tolerance = 1e-10)
  expect_output(print(fit), "survseg")
})

test_that("survseg warns when X column 1 is not the intercept, and coef/summary work", {
  skip_if_not_installed("survival")
  set.seed(4)
  n <- 400; bp <- 200; X <- matrix(1, n, 1)
  d <- sim_survseg(X, par = c(0.5, 3.0, 1.3, 1.3), bp = bp,
                   cens_scale = c(80, 80))
  # First column is not the intercept (constant 2), second is a genuine
  # covariate so the design stays full rank (avoids a survreg rank-deficiency
  # that would abort before the warning path).
  Xbad <- cbind(2, runif(n))
  expect_warning(survseg(d$time, d$status, Xbad, K = 2, init = "binseg"),
                 "intercept")
  fit <- survseg(d$time, d$status, X, K = 2, init = "binseg")
  expect_length(coef(fit), 2L)
  expect_output(summary(fit), "shape")
})
