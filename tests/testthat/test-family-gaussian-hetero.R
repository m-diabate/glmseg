test_that("gaussian_hetero fit matches lm MLE variance and logLik", {
  set.seed(3)
  n <- 300; bp <- 150
  X <- matrix(1, n, 1)
  y <- c(rnorm(bp, 2, 2), rnorm(n - bp, -1, 4))
  fam <- gaussian_hetero()
  fit <- fam$fit_segments(y, X, bp)

  seg1 <- 1:bp; seg2 <- (bp + 1):n
  # MLE sigma per segment: sqrt(RSS_k / n_k)
  s1 <- sqrt(mean((y[seg1] - mean(y[seg1]))^2))
  s2 <- sqrt(mean((y[seg2] - mean(y[seg2]))^2))
  expect_equal(fit$sigma[1], s1, tolerance = 1e-8)
  expect_equal(fit$sigma[2], s2, tolerance = 1e-8)

  # segmentation loglik == sum of per-segment lm logLik
  ll_ref <- as.numeric(logLik(lm(y[seg1] ~ 1))) + as.numeric(logLik(lm(y[seg2] ~ 1)))
  expect_equal(segmentation_loglik(fam, y, X, bp), ll_ref, tolerance = 1e-6)
})

test_that("gaussian_hetero counts (p+1)*K parameters", {
  expect_equal(gaussian_hetero()$n_params(K = 3, p = 2), (2 + 1) * 3)
})
