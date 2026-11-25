test_that("survseg fused-lasso (cox) proposes candidates near a scale change", {
  skip_if_not_installed("survival")
  set.seed(5)
  n <- 500; bp_true <- 300; X <- matrix(1, n, 1)
  d <- sim_survseg(X, par = c(0.3, 3.0, 1.3, 1.3), bp = bp_true,
                   cens_scale = c(80, 80))
  cand <- survseg_init_fused_lasso(d$time, d$status, K = 2)
  expect_true(length(cand) >= 1)
  expect_true(min(abs(cand - bp_true)) <= 70)
})
