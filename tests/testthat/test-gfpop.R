test_that("gfpop_gauss finds a clean single breakpoint", {
  skip_if_not_installed("gfpop")
  set.seed(12)
  n <- 400; bp_true <- 250
  y <- c(rnorm(bp_true, 0, 1), rnorm(n - bp_true, 9, 1))
  res <- gfpop_gauss(y, K = 2)
  expect_true(abs(res$bp[1] - bp_true) <= 3)
})
