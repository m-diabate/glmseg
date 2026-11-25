#' Exhaustive single-breakpoint search (reference method)
#'
#' @param y response vector.
#' @param X n by p design matrix (column 1 = intercept).
#' @param family a [glmseg_family].
#' @param range two quantile probabilities bounding the candidate positions.
#' @return list with `bp` (best breakpoint) and `loglik`.
#' @examples
#' set.seed(1)
#' X <- matrix(1, 400, 1)
#' y <- c(rnorm(250, 0, 1), rnorm(150, 8, 1))
#' bruteforce(y, X, gaussian_homo())$bp
#' @export
bruteforce <- function(y, X, family, range = c(0.05, 0.95)) {
  n <- length(y)
  q <- round(unname(stats::quantile(seq_len(n), range)))
  cand <- q[1L]:q[2L]
  ll <- vapply(cand, function(b) segmentation_loglik(family, y, X, b), numeric(1))
  list(bp = cand[which.max(ll)], loglik = max(ll))
}
