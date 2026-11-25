test_that("survseg_bruteforce finds a clear scale change", {
  skip_if_not_installed("survival")
  set.seed(3)
  n <- 500; bp_true <- 300; X <- matrix(1, n, 1)
  d <- sim_survseg(X, par = c(0.5, 3.0, 1.3, 1.3), bp = bp_true,
                   cens_scale = c(80, 80))
  res <- survseg_bruteforce(d$time, d$status, X)
  expect_true(abs(res$bp - bp_true) <= 40)
})

test_that("survseg_init_binseg returns sorted candidates near breakpoints", {
  skip_if_not_installed("survival")
  set.seed(4)
  n <- 600; X <- matrix(1, n, 1)
  d <- sim_survseg(X, par = c(0.3, 2.5, 0.3, 1.3, 1.3, 1.3),
                   bp = c(200, 400), cens_scale = c(80, 80, 80))
  cand <- survseg_init_binseg(d$time, d$status, X, depth = 3L)
  expect_true(!is.unsorted(cand))
  expect_true(min(abs(cand - 200)) <= 40)
  expect_true(min(abs(cand - 400)) <= 40)
})
