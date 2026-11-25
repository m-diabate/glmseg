make_le <- function(y, mus, sigma) {
  K <- length(mus)
  sapply(seq_len(K), function(k) dnorm(y, mus[k], sigma, log = TRUE))
}

test_that("maxfwbw recovers a clean single breakpoint", {
  set.seed(1)
  n <- 200; bp_true <- 120
  y <- c(rnorm(bp_true, 0, 1), rnorm(n - bp_true, 6, 1))
  le <- make_le(y, c(0, 6), sigma = 1)
  res <- maxfwbw_seg(le, matPijk(2))
  expect_equal(res$bp, bp_true)
  expect_length(res$seg, n)
  expect_equal(res$seg[1], 1L)
  expect_equal(res$seg[n], 2L)
})

test_that("engine criterion equals the hard-segmentation log-likelihood", {
  set.seed(2)
  n <- 200; bp_true <- 120
  y <- c(rnorm(bp_true, 0, 1), rnorm(n - bp_true, 6, 1))
  le <- make_le(y, c(0, 6), sigma = 1)
  res <- maxfwbw_seg(le, matPijk(2))
  hard <- sum(le[cbind(seq_len(n), res$seg)])
  expect_equal(res$loglik, hard, tolerance = 1e-8)
  expect_equal(maxfw_loglik(le, matPijk(2)), res$loglik, tolerance = 1e-8)
})

test_that("maxfwbw enforces R_n = K and stays consistent when the free optimum would not end in K", {
  # State 1 is weakly preferred everywhere; an unconstrained best path would
  # stay in state 1 and never reach K. The constraint R_n = K must still hold,
  # and the reported loglik must equal the emissions summed along the returned
  # segmentation.
  n <- 10
  le <- cbind(rep(0, n), rep(-0.5, n))   # K = 2
  res <- maxfwbw_seg(le, matPijk(2))
  expect_equal(res$seg[n], 2L)
  expect_equal(length(res$bp), 1L)
  hard <- sum(le[cbind(seq_len(n), res$seg)])
  expect_equal(res$loglik, hard, tolerance = 1e-8)
  expect_equal(maxfw_loglik(le, matPijk(2)), res$loglik, tolerance = 1e-8)
})
