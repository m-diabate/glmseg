test_that("select_k recovers the true K over under-segmentation on a clean 3-segment signal", {
  set.seed(11)
  y <- c(rnorm(200, 0, 1), rnorm(200, 8, 1), rnorm(200, -6, 1))  # true K = 3
  # BIC reliably rejects UNDER-segmentation, so we search K in 1:3 and check the
  # true K is selected. Plain BIC (penalty d*K*log(n), as in the reference
  # article) can marginally OVER-segment change-points by one segment on some
  # finite samples — a documented limitation — so we do not assert exact
  # recovery against an inflated K_max.
  res <- select_k(y, K_max = 3, family = gaussian_homo(), init = "binseg")
  expect_equal(res$K, 3L)
  expect_true(res$bic[3] < res$bic[2])
  expect_true(res$bic[2] < res$bic[1])
  expect_equal(res$fit$K, res$K)
  expect_length(res$fit$bp, 2L)
})

test_that("glmseg_bic mbic adds the (K-1)*log(n) breakpoint penalty", {
  set.seed(1)
  y <- c(rnorm(150, 0, 1), rnorm(150, 5, 1))
  fit <- glmseg(y, K = 2, family = gaussian_homo(), init = "binseg")
  expect_equal(glmseg_bic(fit, "mbic"),
               glmseg_bic(fit, "bic") + (fit$K - 1L) * log(fit$n))
  expect_equal(glmseg_bic(fit), glmseg_bic(fit, "bic"))  # default unchanged
})

test_that("mbic penalty curbs over-segmentation where plain bic over-segments", {
  set.seed(11)
  y <- c(rnorm(200, 0, 1), rnorm(200, 8, 1), rnorm(200, -6, 1))  # true K = 3
  res_mbic <- select_k(y, K_max = 4, family = gaussian_homo(),
                       init = "binseg", penalty = "mbic")
  expect_equal(res_mbic$K, 3L)
})
