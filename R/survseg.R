# Fit a Weibull-AFT model per segment. Returns per-segment AFT coefficients
# (first = intercept) and shape (1/scale), or NULL if any segment is unfittable.
survseg_fit_segments <- function(time, status, X, bp) {
  n <- length(time); K <- length(bp) + 1L; bp0 <- c(0L, bp, n)
  betas <- vector("list", K); shape <- numeric(K)
  for (k in seq_len(K)) {
    idx <- (bp0[k] + 1L):bp0[k + 1L]
    Xk <- X[idx, , drop = FALSE]
    df <- data.frame(.time = time[idx], .status = status[idx])
    fit <- suppressWarnings(tryCatch(
      survival::survreg(survival::Surv(.time, .status) ~ Xk - 1,
                        data = df, dist = "weibull"),
      error = function(e) NULL))
    if (is.null(fit) || fit$scale <= 0 || any(!is.finite(fit$coefficients))) {
      return(NULL)
    }
    betas[[k]] <- unname(fit$coefficients)
    shape[k] <- 1 / fit$scale
  }
  list(betas = betas, shape = shape)
}

# n x K matrix of Weibull-AFT log-emissions (Delta*log h + log S).
survseg_log_emission <- function(time, status, X, betas, shape) {
  K <- length(betas); n <- length(time); p <- ncol(X)
  le <- matrix(NA_real_, n, K)
  for (k in seq_len(K)) {
    a <- betas[[k]]; pk <- shape[k]; a0 <- a[1L]
    logh0 <- log(pk) + (pk - 1) * log(time) - pk * a0
    if (p > 1L) {
      betaCov <- -pk * a[-1L]
      XbetaCov <- as.vector(X[, -1L, drop = FALSE] %*% betaCov)
      logh <- logh0 + XbetaCov
      logSurv <- -(time^pk) * exp(-pk * a0) * exp(XbetaCov)
    } else {
      logh <- logh0
      logSurv <- -(time^pk) * exp(-pk * a0)
    }
    le[, k] <- status * logh + logSurv
  }
  le
}

# Log-likelihood of the hard segmentation defined by bp.
survseg_segmentation_loglik <- function(time, status, X, bp) {
  fit <- survseg_fit_segments(time, status, X, bp)
  if (is.null(fit)) return(-Inf)
  le <- survseg_log_emission(time, status, X, fit$betas, fit$shape)
  n <- length(time); bp0 <- c(0L, bp, n)
  seg <- rep.int(seq_along(fit$betas), diff(bp0))
  sum(le[cbind(seq_len(n), seg)])
}

# Exhaustive single-breakpoint search for survival (reference / binseg building block).
survseg_bruteforce <- function(time, status, X, range = c(0.05, 0.95)) {
  n <- length(time)
  q <- round(unname(stats::quantile(seq_len(n), range)))
  cand <- q[1L]:q[2L]
  ll <- vapply(cand,
               function(b) survseg_segmentation_loglik(time, status, X, b),
               numeric(1))
  list(bp = cand[which.max(ll)], loglik = max(ll))
}

# Internal: one max-EM run from a fixed initial breakpoint set (survival).
survseg_maxem <- function(time, status, X, K, lpi, bp, control) {
  n <- length(time); ll_prev <- -Inf
  for (it in seq_len(control$maxit)) {
    fit <- survseg_fit_segments(time, status, X, bp)
    if (is.null(fit)) return(NULL)
    le <- survseg_log_emission(time, status, X, fit$betas, fit$shape)
    if (any(!is.finite(le))) return(NULL)
    step <- maxfwbw_seg(le, lpi)
    if (length(step$bp) != (K - 1L)) return(NULL)
    if (min(diff(c(0L, step$bp, n))) <= ncol(X)) return(NULL)
    if (identical(step$bp, bp) || abs(step$loglik - ll_prev) < control$tol) {
      return(list(bp = step$bp, loglik = step$loglik))
    }
    bp <- step$bp; ll_prev <- step$loglik
  }
  list(bp = bp, loglik = ll_prev)
}

#' Change-point detection for right-censored (Weibull-AFT) data via max-EM
#'
#' @param time vector of observed times.
#' @param status event indicator (1 = event, 0 = censored).
#' @param X optional n by p design matrix (column 1 = intercept). If `NULL`,
#'   an intercept-only matrix is used.
#' @param K number of segments (`K-1` breakpoints).
#' @param init initialization strategy: `"binseg"` or `"fused_lasso"`.
#' @param control see [glmseg_control()].
#' @return an object of class `survseg`.
#' @importFrom gtools combinations
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 300
#'   tt <- c(rweibull(150, 1.3, exp(0.3)), rweibull(150, 1.3, exp(2.5)))
#'   cc <- rexp(n, 1 / 50)
#'   time <- pmin(tt, cc); status <- as.numeric(tt <= cc)
#'   survseg(time, status, K = 2, init = "binseg")$bp
#' }
#' @export
survseg <- function(time, status, X = NULL, K,
                    init = c("binseg", "fused_lasso"),
                    control = glmseg_control()) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for survseg(); please install it.",
         call. = FALSE)
  }
  init <- match.arg(init)
  n <- length(time)
  if (is.null(X)) X <- matrix(1, n, 1L)
  p <- ncol(X)
  if (!all(X[, 1] == 1)) {
    warning("survseg: column 1 of X is expected to be the intercept (all ones); ",
            "results assume this convention.", call. = FALSE)
  }

  if (K == 1L) {
    fit <- survseg_fit_segments(time, status, X, integer(0))
    if (is.null(fit)) stop("survseg: the single-segment model is unfittable.")
    return(structure(
      list(bp = integer(0), betas = fit$betas, shape = fit$shape,
           loglik = survseg_segmentation_loglik(time, status, X, integer(0)),
           K = 1L, n = n, p = p, init = init),
      class = "survseg"))
  }

  lpi <- matPijk(K)
  cand <- if (init == "fused_lasso") {
    survseg_init_fused_lasso(time, status, K, min_seg = control$min_seg)
  } else {
    survseg_init_binseg(time, status, X)
  }
  if (length(cand) < (K - 1L)) {
    stop("Not enough breakpoint candidates for K = ", K,
         " (found ", length(cand), ").")
  }

  combos <- gtools::combinations(length(cand), K - 1L, v = cand,
                                 set = TRUE, repeats.allowed = FALSE)
  best <- NULL
  for (r in seq_len(nrow(combos))) {
    res <- survseg_maxem(time, status, X, K, lpi,
                         as.integer(combos[r, ]), control)
    if (!is.null(res) && (is.null(best) || res$loglik > best$loglik)) {
      best <- res
    }
  }
  if (is.null(best)) stop("All initializations failed to converge.")

  fit <- survseg_fit_segments(time, status, X, best$bp)
  structure(
    list(bp = best$bp, betas = fit$betas, shape = fit$shape,
         loglik = best$loglik, K = K, n = n, p = p, init = init),
    class = "survseg")
}

#' BIC of a fitted survseg model
#'
#' @param fit a `survseg` object.
#' @param penalty `"bic"` (default) or `"mbic"` (adds `(K-1)*log(n)` to curb
#'   over-segmentation).
#' @return the BIC value.
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 300
#'   tt <- c(rweibull(150, 1.3, exp(0.3)), rweibull(150, 1.3, exp(2.5)))
#'   cc <- rexp(n, 1 / 50)
#'   time <- pmin(tt, cc); status <- as.numeric(tt <= cc)
#'   fit <- survseg(time, status, K = 2, init = "binseg")
#'   survseg_bic(fit)
#' }
#' @export
survseg_bic <- function(fit, penalty = c("bic", "mbic")) {
  penalty <- match.arg(penalty)
  np <- (fit$p + 1L) * fit$K
  extra <- if (penalty == "mbic") (fit$K - 1L) else 0L
  -2 * fit$loglik + (np + extra) * log(fit$n)
}

#' Select the number of segments for survseg by BIC
#'
#' @param time vector of observed times.
#' @param status event indicator (1 = event, 0 = censored).
#' @param X optional n by p design matrix (column 1 = intercept).
#' @param K_max maximum number of segments to try.
#' @param init initialization strategy.
#' @param control see [glmseg_control()].
#' @param penalty `"bic"` (default) or `"mbic"`.
#' @return list with the best `fit`, chosen `K`, and the `bic` vector.
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 300
#'   tt <- c(rweibull(150, 1.3, exp(0.3)), rweibull(150, 1.3, exp(2.5)))
#'   cc <- rexp(n, 1 / 50)
#'   time <- pmin(tt, cc); status <- as.numeric(tt <= cc)
#'   survseg_select_k(time, status, K_max = 2, init = "binseg")$K
#' }
#' @export
survseg_select_k <- function(time, status, X = NULL, K_max,
                             init = "binseg", control = glmseg_control(),
                             penalty = c("bic", "mbic")) {
  penalty <- match.arg(penalty)
  fits <- vector("list", K_max); bic <- rep(Inf, K_max)
  for (K in seq_len(K_max)) {
    fit <- tryCatch(survseg(time, status, X, K = K, init = init,
                            control = control), error = function(e) NULL)
    if (!is.null(fit)) { fits[[K]] <- fit; bic[K] <- survseg_bic(fit, penalty) }
  }
  if (all(is.infinite(bic))) {
    stop("survseg_select_k: no value of K produced a valid fit.", call. = FALSE)
  }
  k_best <- which.min(bic)
  list(fit = fits[[k_best]], K = as.integer(k_best), bic = bic)
}
