test_that("poisson_family emission matches dpois", {
  set.seed(1); n <- 60
  X <- cbind(1, rnorm(n)); beta <- c(0.5, 0.3)
  y <- rpois(n, exp(as.vector(X %*% beta)))
  le <- poisson_family()$log_emission(y, X, list(beta), NULL)
  expect_equal(le[, 1], dpois(y, exp(as.vector(X %*% beta)), log = TRUE))
})

test_that("poisson_family counts p*K params and reports poisson glmnet family", {
  expect_equal(poisson_family()$n_params(4, 1), 4)
  expect_equal(poisson_family()$glmnet_family, "poisson")
})

test_that("glmseg recovers a poisson breakpoint with binseg and fused_lasso", {
  set.seed(21); n <- 600; bp_true <- 300; X <- matrix(1, n, 1)
  y <- c(rpois(bp_true, exp(0.4)), rpois(n - bp_true, exp(2.2)))
  fam <- poisson_family()
  fit_bs <- glmseg(y, X, K = 2, family = fam, init = "binseg")
  expect_length(fit_bs$bp, 1L)
  expect_true(abs(fit_bs$bp - bp_true) <= 40)
  fit_fl <- glmseg(y, X, K = 2, family = fam, init = "fused_lasso")
  expect_length(fit_fl$bp, 1L)
  expect_true(abs(fit_fl$bp - bp_true) <= 60)
})

test_that("select_k recovers K over under-segmentation (poisson)", {
  set.seed(22); n <- 900; X <- matrix(1, n, 1)
  y <- c(rpois(300, exp(0.3)), rpois(300, exp(2.3)), rpois(300, exp(0.3)))
  res <- select_k(y, X, K_max = 3, family = poisson_family(), init = "binseg")
  expect_equal(res$K, 3L)
})
