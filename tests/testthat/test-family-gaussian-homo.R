pooled_rss <- function(y, X, bp) {
  n <- length(y); bp0 <- c(0L, bp, n); tot <- 0
  for (k in seq_len(length(bp) + 1L)) {
    idx <- (bp0[k] + 1L):bp0[k + 1L]
    tot <- tot + sum(stats::lm.fit(X[idx, , drop = FALSE], y[idx])$residuals^2)
  }
  tot
}

test_that("gaussian_homo uses a single pooled MLE sigma", {
  set.seed(4)
  n <- 300; bp <- 150; X <- matrix(1, n, 1)
  y <- c(rnorm(bp, 0, 3), rnorm(n - bp, 5, 3))
  fit <- gaussian_homo()$fit_segments(y, X, bp)
  expect_length(fit$sigma, 1L)
  expect_equal(fit$sigma, sqrt(pooled_rss(y, X, bp) / n), tolerance = 1e-10)
})

test_that("unbiased homoscedastic sigma df depends on K (n - K*p), fixing the K>=3 bug", {
  set.seed(5)
  n <- 600; bp <- c(200, 400); X <- matrix(1, n, 1)  # K = 3, p = 1
  y <- c(rnorm(200, 0, 2), rnorm(200, 4, 2), rnorm(200, -3, 2))
  fit <- gaussian_homo()$fit_segments(y, X, bp)
  rss <- pooled_rss(y, X, bp); K <- 3; p <- 1
  expect_equal(fit$sigma_unbiased, sqrt(rss / (n - K * p)), tolerance = 1e-10)
  # the old buggy formula divided by n-(p+1)=n-2, which differs for K>=3:
  expect_false(isTRUE(all.equal(fit$sigma_unbiased, sqrt(rss / (n - (p + 1))))))
})

test_that("gaussian_homo counts K*p + 1 parameters", {
  expect_equal(gaussian_homo()$n_params(K = 3, p = 2), 3 * 2 + 1)
})
