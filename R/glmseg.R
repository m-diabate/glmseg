# Internal: in-sample fitted values (per-segment linear predictor by position).
glmseg_fitted_vals <- function(X, betas, bp) {
  n <- nrow(X)
  bp0 <- c(0L, bp, n)
  out <- numeric(n)
  for (k in seq_along(betas)) {
    idx <- (bp0[k] + 1L):bp0[k + 1L]
    out[idx] <- as.vector(X[idx, , drop = FALSE] %*% betas[[k]])
  }
  out
}

#' Control parameters for [glmseg()]
#'
#' @param maxit maximum max-EM iterations per initialization.
#' @param tol convergence tolerance on the criterion.
#' @param min_seg minimum individuals per segment (used by fused-lasso).
#' @param verbose emit progress messages.
#' @return a list of control parameters.
#' @examples
#' ctrl <- glmseg_control(maxit = 50, tol = 1e-6)
#' set.seed(1)
#' y <- c(rnorm(150, 0, 1), rnorm(150, 6, 1))
#' glmseg(y, K = 2, family = gaussian_homo(), init = "binseg", control = ctrl)$bp
#' @export
glmseg_control <- function(maxit = 100L, tol = 1e-5, min_seg = 50L,
                           verbose = FALSE) {
  list(maxit = maxit, tol = tol, min_seg = min_seg, verbose = verbose)
}

#' Change-point detection via the max-EM algorithm
#'
#' @param y response vector.
#' @param X optional n by p design matrix (column 1 = intercept). If `NULL`, an
#'   intercept-only matrix is used.
#' @param K number of segments (`K-1` breakpoints).
#' @param family a [glmseg_family], e.g. [gaussian_hetero()] or [gaussian_homo()].
#' @param init initialization strategy: `"fused_lasso"` or `"binseg"`.
#' @param control see [glmseg_control()].
#' @return an object of class `glmseg`.
#' @importFrom gtools combinations
#' @examples
#' set.seed(1)
#' y <- c(rnorm(200, 0, 1), rnorm(200, 8, 1), rnorm(200, -6, 1))
#' fit <- glmseg(y, K = 3, family = gaussian_homo(), init = "binseg")
#' fit$bp
#' summary(fit)
#' @export
glmseg <- function(y, X = NULL, K, family = gaussian_hetero(),
                   init = c("fused_lasso", "binseg"),
                   control = glmseg_control()) {
  init <- match.arg(init)
  n <- length(y)
  if (is.null(X)) X <- matrix(1, n, 1L)
  p <- ncol(X)
  if (!all(X[, 1] == 1)) {
    warning("glmseg: column 1 of X is expected to be the intercept (all ones); ",
            "results assume this convention.", call. = FALSE)
  }
  if (K == 1L) {
    fit <- family$fit_segments(y, X, integer(0))
    fitted_vals <- glmseg_fitted_vals(X, fit$betas, integer(0))
    return(structure(
      list(bp = integer(0), betas = fit$betas, sigma = fit$sigma,
           sigma_unbiased = fit$sigma_unbiased,
           fitted = fitted_vals,
           loglik = segmentation_loglik(family, y, X, integer(0)),
           K = 1L, family = family, init = init, n = n, p = p),
      class = "glmseg"))
  }
  lpi <- matPijk(K)

  cand <- if (init == "fused_lasso") {
    init_fused_lasso(y, K, glmnet_family = family$glmnet_family,
                     min_seg = control$min_seg)
  } else {
    init_binseg(y, X, family)
  }
  if (length(cand) < (K - 1L)) {
    stop("Not enough breakpoint candidates for K = ", K,
         " (found ", length(cand), ").")
  }

  combos <- gtools::combinations(length(cand), K - 1L, v = cand,
                                 set = TRUE, repeats.allowed = FALSE)
  best <- NULL
  for (r in seq_len(nrow(combos))) {
    bp <- as.integer(combos[r, ])
    res <- maxem_iterate(y, X, K, family, lpi, bp, control)
    if (!is.null(res) && (is.null(best) || res$loglik > best$loglik)) {
      best <- res
    }
  }
  if (is.null(best)) stop("All initializations failed to converge.")

  fit <- family$fit_segments(y, X, best$bp)
  fitted_vals <- glmseg_fitted_vals(X, fit$betas, best$bp)
  structure(
    list(bp = best$bp, betas = fit$betas, sigma = fit$sigma,
         sigma_unbiased = fit$sigma_unbiased, fitted = fitted_vals,
         loglik = best$loglik,
         K = K, family = family, init = init, n = n, p = p),
    class = "glmseg"
  )
}

# Internal: run max-EM from a fixed initial breakpoint set.
maxem_iterate <- function(y, X, K, family, lpi, bp, control) {
  n <- length(y)
  ll_prev <- -Inf
  for (it in seq_len(control$maxit)) {
    fit <- family$fit_segments(y, X, bp)
    le <- family$log_emission(y, X, fit$betas, fit$sigma)
    step <- maxfwbw_seg(le, lpi)
    if (length(step$bp) != (K - 1L)) return(NULL)          # degenerate
    if (min(diff(c(0L, step$bp, n))) <= ncol(X)) return(NULL)   # need > p points per segment
    if (identical(step$bp, bp) || abs(step$loglik - ll_prev) < control$tol) {
      return(list(bp = step$bp, loglik = step$loglik))
    }
    bp <- step$bp
    ll_prev <- step$loglik
  }
  list(bp = bp, loglik = ll_prev)
}
