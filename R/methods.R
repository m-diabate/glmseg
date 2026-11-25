#' Print a glmseg fit
#'
#' @param x a `glmseg` object.
#' @param ... ignored.
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' y <- c(rnorm(150, 0, 1), rnorm(150, 5, 1))
#' fit <- glmseg(y, K = 2, family = gaussian_homo(), init = "binseg")
#' print(fit)
#' @export
print.glmseg <- function(x, ...) {
  cat("glmseg fit (", x$family$name, ")\n", sep = "")
  cat("  K = ", x$K, " segments, breakpoints: ",
      if (length(x$bp)) paste(x$bp, collapse = ", ") else "none", "\n", sep = "")
  cat("  logLik = ", format(x$loglik), ",  BIC = ", format(glmseg_bic(x)), "\n",
      sep = "")
  invisible(x)
}

#' Summarize a glmseg fit
#'
#' @param object a `glmseg` object.
#' @param ... ignored.
#' @return `object`, invisibly.
#' @examples
#' set.seed(1)
#' y <- c(rnorm(150, 0, 1), rnorm(150, 5, 1))
#' fit <- glmseg(y, K = 2, family = gaussian_homo(), init = "binseg")
#' summary(fit)
#' @export
summary.glmseg <- function(object, ...) {
  cat("glmseg fit (", object$family$name, ")\n", sep = "")
  cat("  breakpoints:", if (length(object$bp)) paste(object$bp, collapse = ", ")
      else "none", "\n")
  for (k in seq_along(object$betas)) {
    cat("  segment", k, "beta:", paste(format(object$betas[[k]]), collapse = " "),
        "\n")
  }
  if (!is.null(object$sigma)) {
    cat("  sigma (MLE):", paste(format(object$sigma), collapse = " "), "\n")
    cat("  sigma (unbiased):",
        paste(format(object$sigma_unbiased), collapse = " "), "\n")
  }
  invisible(object)
}

#' Segment regression coefficients of a glmseg fit
#'
#' @param object a `glmseg` object.
#' @param ... ignored.
#' @return a list of per-segment coefficient vectors.
#' @importFrom stats coef
#' @examples
#' set.seed(1)
#' y <- c(rnorm(150, 0, 1), rnorm(150, 5, 1))
#' fit <- glmseg(y, K = 2, family = gaussian_homo(), init = "binseg")
#' coef(fit)
#' @export
coef.glmseg <- function(object, ...) object$betas

#' Log-likelihood of a glmseg fit
#'
#' @param object a `glmseg` object.
#' @param ... ignored.
#' @return a `logLik` object with `df` and `nobs` attributes.
#' @importFrom stats logLik
#' @examples
#' set.seed(1)
#' y <- c(rnorm(150, 0, 1), rnorm(150, 5, 1))
#' fit <- glmseg(y, K = 2, family = gaussian_homo(), init = "binseg")
#' logLik(fit)
#' @export
logLik.glmseg <- function(object, ...) {
  val <- object$loglik
  attr(val, "df") <- object$family$n_params(object$K, object$p)
  attr(val, "nobs") <- object$n
  class(val) <- "logLik"
  val
}

#' BIC of a glmseg fit
#'
#' @param object a `glmseg` object.
#' @param ... ignored.
#' @return the BIC value.
#' @importFrom stats BIC
#' @examples
#' set.seed(1)
#' y <- c(rnorm(150, 0, 1), rnorm(150, 5, 1))
#' fit <- glmseg(y, K = 2, family = gaussian_homo(), init = "binseg")
#' BIC(fit)
#' @export
BIC.glmseg <- function(object, ...) glmseg_bic(object)

#' Fitted values of a glmseg fit
#'
#' @param object a `glmseg` object.
#' @param ... ignored.
#' @return the in-sample fitted values (numeric vector of length `n`).
#' @importFrom stats fitted
#' @examples
#' set.seed(1)
#' y <- c(rnorm(150, 0, 1), rnorm(150, 5, 1))
#' fit <- glmseg(y, K = 2, family = gaussian_homo(), init = "binseg")
#' fitted(fit)
#' @export
fitted.glmseg <- function(object, ...) {
  object$fitted
}

#' Predictions from a glmseg fit
#'
#' @param object a `glmseg` object.
#' @param X optional design matrix aligned by position with the fitted
#'   segmentation (same rows/order). If `NULL`, the in-sample fitted values
#'   are returned.
#' @param ... ignored.
#' @return a numeric vector of per-segment predictions.
#' @importFrom stats predict
#' @examples
#' set.seed(1)
#' y <- c(rnorm(150, 0, 1), rnorm(150, 5, 1))
#' fit <- glmseg(y, K = 2, family = gaussian_homo(), init = "binseg")
#' predict(fit)
#' @export
predict.glmseg <- function(object, X = NULL, ...) {
  if (is.null(X)) return(object$fitted)
  glmseg_fitted_vals(X, object$betas, object$bp)
}

#' Plot a glmseg fit
#'
#' @param x a `glmseg` object.
#' @param y the observed response vector to plot.
#' @param ... passed to the underlying plot call.
#' @return `x`, invisibly.
#' @importFrom graphics plot
#' @examples
#' set.seed(1)
#' y <- c(rnorm(150, 0, 1), rnorm(150, 5, 1))
#' fit <- glmseg(y, K = 2, family = gaussian_homo(), init = "binseg")
#' plot(fit, y = y)
#' @export
plot.glmseg <- function(x, y = NULL, ...) {
  if (is.null(y)) stop("Supply the observed response via 'y' to plot.")
  plot(y, type = "p", pch = 20, col = "grey50",
       xlab = "index", ylab = "y", ...)
  graphics::lines(fitted(x), col = "red", lwd = 2)
  if (length(x$bp)) graphics::abline(v = x$bp, lty = 2, col = "blue")
  invisible(x)
}
