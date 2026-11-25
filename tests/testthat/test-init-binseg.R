test_that("binary segmentation returns sorted candidates covering breakpoints", {
  set.seed(8)
  n <- 600; X <- matrix(1, n, 1)
  y <- c(rnorm(200, 0, 1), rnorm(200, 6, 1), rnorm(200, -4, 1))
  cand <- init_binseg(y, X, gaussian_homo(), depth = 3L)
  expect_true(!is.unsorted(cand))
  expect_true(min(abs(cand - 200)) <= 20)
  expect_true(min(abs(cand - 400)) <= 20)
})
