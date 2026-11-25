test_that("logistic_family emission matches dbinom", {
  set.seed(1); n <- 60
  X <- cbind(1, rnorm(n)); beta <- c(0.3, -0.8)
  y <- rbinom(n, 1, plogis(as.vector(X %*% beta)))
  le <- logistic_family()$log_emission(y, X, list(beta), NULL)
  expect_equal(le[, 1],
               dbinom(y, 1, plogis(as.vector(X %*% beta)), log = TRUE))
})

test_that("logistic_family counts p*K params and reports binomial glmnet family", {
  expect_equal(logistic_family()$n_params(3, 2), 6)
  expect_equal(logistic_family()$glmnet_family, "binomial")
})

test_that("glmseg recovers a logistic breakpoint with binseg and fused_lasso", {
  set.seed(11); n <- 600; bp_true <- 300; X <- matrix(1, n, 1)
  y <- c(rbinom(bp_true, 1, plogis(-1.6)),
         rbinom(n - bp_true, 1, plogis(1.6)))
  fam <- logistic_family()
  fit_bs <- glmseg(y, X, K = 2, family = fam, init = "binseg")
  expect_length(fit_bs$bp, 1L)
  expect_true(abs(fit_bs$bp - bp_true) <= 60)
  fit_fl <- glmseg(y, X, K = 2, family = fam, init = "fused_lasso")
  expect_length(fit_fl$bp, 1L)
  expect_true(abs(fit_fl$bp - bp_true) <= 80)
})

test_that("select_k recovers K over under-segmentation (logistic)", {
  set.seed(12); n <- 900; X <- matrix(1, n, 1)
  y <- c(rbinom(300, 1, plogis(-1.8)),
         rbinom(300, 1, plogis( 1.8)),
         rbinom(300, 1, plogis(-1.8)))
  res <- select_k(y, X, K_max = 3, family = logistic_family(), init = "binseg")
  expect_equal(res$K, 3L)
})
