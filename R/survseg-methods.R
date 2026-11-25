#' Print a survseg fit
#' @param x a `survseg` object.
#' @param ... ignored.
#' @return `x`, invisibly.
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 300
#'   tt <- c(rweibull(150, 1.3, exp(0.3)), rweibull(150, 1.3, exp(2.5)))
#'   cc <- rexp(n, 1 / 50)
#'   time <- pmin(tt, cc); status <- as.numeric(tt <= cc)
#'   fit <- survseg(time, status, K = 2, init = "binseg")
#'   print(fit)
#' }
#' @export
print.survseg <- function(x, ...) {
  cat("survseg fit (Weibull-AFT)\n")
  cat("  K = ", x$K, " segments, breakpoints: ",
      if (length(x$bp)) paste(x$bp, collapse = ", ") else "none", "\n", sep = "")
  cat("  logLik = ", format(x$loglik), ",  BIC = ", format(survseg_bic(x)),
      "\n", sep = "")
  invisible(x)
}

#' Summarize a survseg fit
#' @param object a `survseg` object.
#' @param ... ignored.
#' @return `object`, invisibly.
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 300
#'   tt <- c(rweibull(150, 1.3, exp(0.3)), rweibull(150, 1.3, exp(2.5)))
#'   cc <- rexp(n, 1 / 50)
#'   time <- pmin(tt, cc); status <- as.numeric(tt <= cc)
#'   fit <- survseg(time, status, K = 2, init = "binseg")
#'   summary(fit)
#' }
#' @export
summary.survseg <- function(object, ...) {
  cat("survseg fit (Weibull-AFT)\n")
  cat("  breakpoints:", if (length(object$bp)) paste(object$bp, collapse = ", ")
      else "none", "\n")
  for (k in seq_along(object$betas)) {
    cat("  segment", k, "AFT beta:",
        paste(format(object$betas[[k]]), collapse = " "),
        "| shape:", format(object$shape[k]), "\n")
  }
  invisible(object)
}

#' Segment AFT coefficients of a survseg fit
#' @param object a `survseg` object.
#' @param ... ignored.
#' @return a list of per-segment AFT coefficient vectors.
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 300
#'   tt <- c(rweibull(150, 1.3, exp(0.3)), rweibull(150, 1.3, exp(2.5)))
#'   cc <- rexp(n, 1 / 50)
#'   time <- pmin(tt, cc); status <- as.numeric(tt <= cc)
#'   fit <- survseg(time, status, K = 2, init = "binseg")
#'   coef(fit)
#' }
#' @export
coef.survseg <- function(object, ...) object$betas

#' Log-likelihood of a survseg fit
#' @param object a `survseg` object.
#' @param ... ignored.
#' @return a `logLik` object with `df` and `nobs` attributes.
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 300
#'   tt <- c(rweibull(150, 1.3, exp(0.3)), rweibull(150, 1.3, exp(2.5)))
#'   cc <- rexp(n, 1 / 50)
#'   time <- pmin(tt, cc); status <- as.numeric(tt <= cc)
#'   fit <- survseg(time, status, K = 2, init = "binseg")
#'   logLik(fit)
#' }
#' @export
logLik.survseg <- function(object, ...) {
  val <- object$loglik
  attr(val, "df") <- (object$p + 1L) * object$K
  attr(val, "nobs") <- object$n
  class(val) <- "logLik"
  val
}

#' BIC of a survseg fit
#' @param object a `survseg` object.
#' @param ... ignored.
#' @return the BIC value.
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 300
#'   tt <- c(rweibull(150, 1.3, exp(0.3)), rweibull(150, 1.3, exp(2.5)))
#'   cc <- rexp(n, 1 / 50)
#'   time <- pmin(tt, cc); status <- as.numeric(tt <= cc)
#'   fit <- survseg(time, status, K = 2, init = "binseg")
#'   BIC(fit)
#' }
#' @export
BIC.survseg <- function(object, ...) survseg_bic(object)

#' Convert survseg AFT coefficients to proportional-hazards (PH) coefficients
#'
#' For a Weibull model, the PH coefficients relate to the AFT coefficients by
#' `beta_PH = -shape * alpha_AFT` (see the AFT/PH note).
#' @param fit a `survseg` object.
#' @return a list of per-segment PH coefficient vectors.
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 300
#'   tt <- c(rweibull(150, 1.3, exp(0.3)), rweibull(150, 1.3, exp(2.5)))
#'   cc <- rexp(n, 1 / 50)
#'   time <- pmin(tt, cc); status <- as.numeric(tt <= cc)
#'   fit <- survseg(time, status, K = 2, init = "binseg")
#'   aft_to_ph(fit)
#' }
#' @export
aft_to_ph <- function(fit) {
  Map(function(a, s) -s * a, fit$betas, fit$shape)
}
