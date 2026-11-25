#' Gaussian family with segment-specific variance (heteroscedastic)
#'
#' @return a [glmseg_family].
#' @examples
#' set.seed(1)
#' y <- c(rnorm(150, 0, 1), rnorm(150, 5, 3))
#' glmseg(y, K = 2, family = gaussian_hetero(), init = "binseg")$bp
#' @export
gaussian_hetero <- function() {
  new_glmseg_family(
    name = "gaussian_hetero",
    fit_segments = function(y, X, bp) {
      n <- length(y); K <- length(bp) + 1L
      bp0 <- c(0L, bp, n)
      betas <- vector("list", K)
      sigma <- numeric(K); sigma_unb <- numeric(K)
      p <- ncol(X)
      for (k in seq_len(K)) {
        idx <- (bp0[k] + 1L):bp0[k + 1L]
        fit <- stats::lm.fit(X[idx, , drop = FALSE], y[idx])
        betas[[k]] <- fit$coefficients
        rss <- sum(fit$residuals^2); nk <- length(idx)
        sigma[k] <- sqrt(rss / nk)               # MLE
        sigma_unb[k] <- sqrt(rss / max(nk - p, 1L))  # unbiased attribute
      }
      list(betas = betas, sigma = sigma, sigma_unbiased = sigma_unb)
    },
    log_emission = function(y, X, betas, sigma) {
      K <- length(betas); n <- length(y)
      le <- matrix(NA_real_, n, K)
      for (k in seq_len(K)) {
        mu <- as.vector(X %*% betas[[k]])
        le[, k] <- stats::dnorm(y, mu, sigma[k], log = TRUE)
      }
      le
    },
    n_params = function(K, p) (p + 1L) * K
  )
}

#' Gaussian family with common variance (homoscedastic)
#'
#' @return a [glmseg_family]. The likelihood uses the pooled MLE variance
#'   `sum_k RSS_k / n`; the unbiased estimator `sum_k RSS_k / (n - K*p)` is
#'   returned as `sigma_unbiased`.
#' @examples
#' set.seed(1)
#' y <- c(rnorm(150, 0, 1), rnorm(150, 5, 1))
#' glmseg(y, K = 2, family = gaussian_homo(), init = "binseg")$bp
#' @export
gaussian_homo <- function() {
  new_glmseg_family(
    name = "gaussian_homo",
    fit_segments = function(y, X, bp) {
      n <- length(y); K <- length(bp) + 1L; p <- ncol(X)
      bp0 <- c(0L, bp, n)
      betas <- vector("list", K); rss <- 0
      for (k in seq_len(K)) {
        idx <- (bp0[k] + 1L):bp0[k + 1L]
        fit <- stats::lm.fit(X[idx, , drop = FALSE], y[idx])
        betas[[k]] <- fit$coefficients
        rss <- rss + sum(fit$residuals^2)
      }
      list(
        betas = betas,
        sigma = sqrt(rss / n),                       # pooled MLE
        sigma_unbiased = sqrt(rss / max(n - K * p, 1L))
      )
    },
    log_emission = function(y, X, betas, sigma) {
      K <- length(betas); n <- length(y)
      le <- matrix(NA_real_, n, K)
      for (k in seq_len(K)) {
        mu <- as.vector(X %*% betas[[k]])
        le[, k] <- stats::dnorm(y, mu, sigma, log = TRUE)  # single shared sigma
      }
      le
    },
    n_params = function(K, p) K * p + 1L
  )
}
