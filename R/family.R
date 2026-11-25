#' Construct a glmseg family object
#'
#' @param name character identifier.
#' @param fit_segments function `(y, X, bp)` returning at least `betas`
#'   (list of length K), optionally `sigma`/`sigma_unbiased`.
#' @param log_emission function `(y, X, betas, sigma)` returning an n by K matrix.
#' @param n_params function `(K, p)` returning the number of estimated
#'   parameters (dispersion included) for BIC.
#' @param glmnet_family the `glmnet` family string used by the Fused-Lasso
#'   initialization (`"gaussian"`, `"binomial"`, `"poisson"`).
#' @return an object of class `glmseg_family`.
#' @aliases glmseg_family
#' @examples
#' # Built-in families are glmseg_family objects:
#' fam <- gaussian_homo()
#' class(fam)
#' fam$glmnet_family
#' @export
new_glmseg_family <- function(name, fit_segments, log_emission, n_params,
                              glmnet_family = "gaussian") {
  stopifnot(is.character(name), is.function(fit_segments),
            is.function(log_emission), is.function(n_params),
            is.character(glmnet_family))
  structure(
    list(name = name, fit_segments = fit_segments,
         log_emission = log_emission, n_params = n_params,
         glmnet_family = glmnet_family),
    class = "glmseg_family"
  )
}

#' Log-likelihood of a hard segmentation
#'
#' @param family a [glmseg_family].
#' @param y response vector.
#' @param X n by p design matrix (column 1 = intercept).
#' @param bp integer breakpoints (may be empty for K = 1).
#' @return scalar log-likelihood of the segmentation defined by `bp`.
#' @examples
#' set.seed(1)
#' X <- matrix(1, 300, 1)
#' y <- c(rnorm(150, 0, 1), rnorm(150, 5, 1))
#' segmentation_loglik(gaussian_homo(), y, X, bp = 150)
#' @export
segmentation_loglik <- function(family, y, X, bp) {
  fit <- family$fit_segments(y, X, bp)
  le <- family$log_emission(y, X, fit$betas, fit$sigma)
  n <- length(y)
  bp0 <- c(0L, bp, n)
  seg <- rep.int(seq_along(fit$betas), diff(bp0))
  sum(le[cbind(seq_len(n), seg)])
}
