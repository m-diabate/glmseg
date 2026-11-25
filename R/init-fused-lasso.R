#' Fused-Lasso breakpoint candidates
#'
#' @param y response vector.
#' @param K target number of segments.
#' @param glmnet_family the `glmnet` family for the penalized fit
#'   (`"gaussian"`, `"binomial"`, `"poisson"`).
#' @param min_seg minimum number of individuals per segment.
#' @param factor_lo multiplier for the lower bound on candidate count.
#' @return integer vector of candidate breakpoints.
#' @importFrom glmnet glmnet
#' @examples
#' set.seed(1)
#' y <- c(rnorm(250, 0, 1), rnorm(250, 10, 1))
#' init_fused_lasso(y, K = 2)
#' @export
init_fused_lasso <- function(y, K, glmnet_family = "gaussian",
                             min_seg = 50L, factor_lo = 5L) {
  n <- length(y)
  D <- matrix(1, n, n); D[upper.tri(D)] <- 0
  fit <- glmnet::glmnet(x = D, y = y, family = glmnet_family)
  n_seg <- apply(fit$beta, 2L, function(b) sum(abs(b) > 0))
  reached <- which(n_seg >= factor_lo * (K - 1L))
  if (length(reached) == 0L) {
    stop("Fused-Lasso initialization did not reach the requested resolution; ",
         "try init = 'binseg' or a different K.", call. = FALSE)
  }
  col <- min(reached)
  cand <- which(abs(fit$beta[, col]) > 0)

  it <- 1L; upper <- (3 * (K - 1L)) / 2
  while (length(cand) > upper && it < min_seg) {
    cand <- almost_unique_tol(cand, tol = it)
    it <- it + 1L
  }
  cand <- cand[cand > 2L & cand < (n - 1L)]
  as.integer(cand)
}
