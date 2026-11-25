test_that("families carry a glmnet_family, defaulting to gaussian", {
  expect_equal(gaussian_homo()$glmnet_family, "gaussian")
  expect_equal(gaussian_hetero()$glmnet_family, "gaussian")
  fam <- new_glmseg_family("dummy",
                           fit_segments = function(y, X, bp) list(betas = list()),
                           log_emission = function(y, X, betas, sigma) matrix(0),
                           n_params = function(K, p) 0L)
  expect_equal(fam$glmnet_family, "gaussian")
  fam2 <- new_glmseg_family("d2",
                            fit_segments = function(y, X, bp) list(betas = list()),
                            log_emission = function(y, X, betas, sigma) matrix(0),
                            n_params = function(K, p) 0L,
                            glmnet_family = "poisson")
  expect_equal(fam2$glmnet_family, "poisson")
})
