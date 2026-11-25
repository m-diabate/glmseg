#' BIC of a fitted glmseg model
#'
#' @param fit a `glmseg` object.
#' @param penalty `"bic"` (default) uses `d*K` parameters (the article
#'   criterion); `"mbic"` additionally penalizes the `K-1` breakpoint locations
#'   (`+(K-1)*log(n)`), which curbs over-segmentation.
#' @return the BIC value.
#' @examples
#' set.seed(1)
#' y <- c(rnorm(150, 0, 1), rnorm(150, 5, 1))
#' fit <- glmseg(y, K = 2, family = gaussian_homo(), init = "binseg")
#' glmseg_bic(fit)
#' glmseg_bic(fit, penalty = "mbic")
#' @export
glmseg_bic <- function(fit, penalty = c("bic", "mbic")) {
  penalty <- match.arg(penalty)
  np <- fit$family$n_params(fit$K, fit$p)
  extra <- if (penalty == "mbic") (fit$K - 1L) else 0L
  -2 * fit$loglik + (np + extra) * log(fit$n)
}

#' Select the number of segments by BIC
#'
#' @param y response vector.
#' @param X optional n by p design matrix (column 1 = intercept).
#' @param K_max maximum number of segments to try.
#' @param family a [glmseg_family].
#' @param init initialization strategy.
#' @param control see [glmseg_control()].
#' @param penalty criterion passed to [glmseg_bic()]: `"bic"` (default) or
#'   `"mbic"` (penalizes breakpoints to curb over-segmentation).
#' @return list with the best `fit`, chosen `K`, and the `bic` vector.
#' @examples
#' set.seed(1)
#' y <- c(rnorm(200, 0, 1), rnorm(200, 8, 1), rnorm(200, -6, 1))
#' res <- select_k(y, K_max = 4, family = gaussian_homo(),
#'                 init = "binseg", penalty = "mbic")
#' res$K
#' @export
select_k <- function(y, X = NULL, K_max, family = gaussian_hetero(),
                     init = "binseg", control = glmseg_control(),
                     penalty = c("bic", "mbic")) {
  penalty <- match.arg(penalty)
  fits <- vector("list", K_max)
  bic <- rep(Inf, K_max)
  for (K in seq_len(K_max)) {
    fit <- tryCatch(glmseg(y, X, K = K, family = family, init = init,
                           control = control), error = function(e) NULL)
    if (!is.null(fit)) { fits[[K]] <- fit; bic[K] <- glmseg_bic(fit, penalty) }
  }
  if (all(is.infinite(bic))) {
    stop("select_k: no value of K produced a valid fit.", call. = FALSE)
  }
  k_best <- which.min(bic)
  list(fit = fits[[k_best]], K = as.integer(k_best), bic = bic)
}
