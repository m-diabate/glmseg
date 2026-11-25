test_that("glmseg S3 accessors return coherent values", {
  set.seed(13)
  n <- 400; X <- matrix(1, n, 1)
  y <- c(rnorm(200, 0, 1), rnorm(200, 6, 1))
  fit <- glmseg(y, X, K = 2, family = gaussian_homo(), init = "binseg")

  expect_s3_class(fit, "glmseg")
  expect_equal(as.numeric(logLik(fit)), fit$loglik)
  expect_equal(attr(logLik(fit), "df"), fit$family$n_params(fit$K, fit$p))
  expect_equal(BIC(fit), -2 * fit$loglik +
                 fit$family$n_params(fit$K, fit$p) * log(n))
  expect_length(fitted(fit), n)
  expect_output(print(fit), "glmseg")
})

test_that("glmseg warns when X column 1 is not the intercept", {
  set.seed(2)
  n <- 300; y <- c(rnorm(150, 0, 1), rnorm(150, 5, 1))
  Xbad <- cbind(2, rnorm(n))          # first column is 2, not 1
  expect_warning(glmseg(y, Xbad, K = 2, family = gaussian_hetero(),
                        init = "binseg"), "intercept")
  Xok <- matrix(1, n, 1)
  expect_silent(glmseg(y, Xok, K = 2, family = gaussian_homo(),
                       init = "binseg"))
})

test_that("coef and summary work for a glmseg fit", {
  set.seed(3)
  n <- 300; y <- c(rnorm(150, 0, 1), rnorm(150, 5, 1))
  fit <- glmseg(y, K = 2, family = gaussian_homo(), init = "binseg")
  cf <- coef(fit)
  expect_length(cf, 2L)
  expect_output(summary(fit), "sigma")
})
