#' Logistic (binomial) loss family for glmseg
#'
#' Segment-wise logistic regression; no dispersion parameter.
#' @return a [glmseg_family].
#' @examples
#' set.seed(1)
#' y <- c(rbinom(200, 1, plogis(-1.5)), rbinom(200, 1, plogis(1.5)))
#' glmseg(y, K = 2, family = logistic_family(), init = "binseg")$bp
#' @export
logistic_family <- function() {
  new_glmseg_family(
    name = "logistic",
    fit_segments = function(y, X, bp) {
      n <- length(y); K <- length(bp) + 1L; bp0 <- c(0L, bp, n)
      betas <- vector("list", K)
      for (k in seq_len(K)) {
        idx <- (bp0[k] + 1L):bp0[k + 1L]
        # Separation warnings are expected/benign during the segmentation
        # search; suppress them to keep output clean.
        fit <- suppressWarnings(
          stats::glm.fit(X[idx, , drop = FALSE], y[idx],
                         family = stats::binomial()))
        betas[[k]] <- fit$coefficients
      }
      list(betas = betas)
    },
    log_emission = function(y, X, betas, sigma) {
      K <- length(betas); n <- length(y)
      le <- matrix(NA_real_, n, K)
      for (k in seq_len(K)) {
        eta <- as.vector(X %*% betas[[k]])
        le[, k] <- stats::dbinom(y, size = 1, prob = stats::plogis(eta),
                                 log = TRUE)
      }
      le
    },
    n_params = function(K, p) p * K,
    glmnet_family = "binomial"
  )
}

#' Poisson loss family for glmseg
#'
#' Segment-wise Poisson (log-link) regression; no dispersion parameter.
#' @return a [glmseg_family].
#' @examples
#' set.seed(1)
#' y <- c(rpois(200, exp(0.4)), rpois(200, exp(2.2)))
#' glmseg(y, K = 2, family = poisson_family(), init = "binseg")$bp
#' @export
poisson_family <- function() {
  new_glmseg_family(
    name = "poisson",
    fit_segments = function(y, X, bp) {
      n <- length(y); K <- length(bp) + 1L; bp0 <- c(0L, bp, n)
      betas <- vector("list", K)
      for (k in seq_len(K)) {
        idx <- (bp0[k] + 1L):bp0[k + 1L]
        fit <- suppressWarnings(
          stats::glm.fit(X[idx, , drop = FALSE], y[idx],
                         family = stats::poisson()))
        betas[[k]] <- fit$coefficients
      }
      list(betas = betas)
    },
    log_emission = function(y, X, betas, sigma) {
      K <- length(betas); n <- length(y)
      le <- matrix(NA_real_, n, K)
      for (k in seq_len(K)) {
        lambda <- exp(as.vector(X %*% betas[[k]]))
        le[, k] <- stats::dpois(y, lambda, log = TRUE)
      }
      le
    },
    n_params = function(K, p) p * K,
    glmnet_family = "poisson"
  )
}
