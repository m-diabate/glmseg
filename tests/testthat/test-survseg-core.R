test_that("survseg emission matches survreg logLik at a fixed segmentation", {
  skip_if_not_installed("survival")
  set.seed(1)
  n <- 400; bp <- 200; X <- matrix(1, n, 1)
  d <- sim_survseg(X, par = c(1.5, 3.0, 1.2, 1.2), bp = bp,
                   cens_scale = c(60, 60))
  ll_emis <- survseg_segmentation_loglik(d$time, d$status, X, bp)

  # independent reference: sum of per-segment survreg logLik
  ll_ref <- 0
  for (rng in list(1:bp, (bp + 1):n)) {
    fit <- survival::survreg(
      survival::Surv(d$time[rng], d$status[rng]) ~ 1, dist = "weibull")
    ll_ref <- ll_ref + as.numeric(logLik(fit))
  }
  expect_equal(ll_emis, ll_ref, tolerance = 1e-4)
})

test_that("survseg_fit_segments returns per-segment betas and shape", {
  skip_if_not_installed("survival")
  set.seed(2)
  n <- 400; bp <- 200; X <- matrix(1, n, 1)
  d <- sim_survseg(X, par = c(1.5, 3.0, 1.2, 1.2), bp = bp,
                   cens_scale = c(60, 60))
  fit <- survseg_fit_segments(d$time, d$status, X, bp)
  expect_length(fit$betas, 2L)
  expect_length(fit$shape, 2L)
  expect_true(all(fit$shape > 0))
})
