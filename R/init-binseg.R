#' Binary Segmentation breakpoint candidates
#'
#' @param y response vector.
#' @param X n by p design matrix (column 1 = intercept).
#' @param family a [glmseg_family].
#' @param depth recursion depth (levels of splitting).
#' @param min_seg minimum segment length to attempt a split.
#' @return sorted integer vector of candidate breakpoints.
#' @examples
#' set.seed(1)
#' X <- matrix(1, 400, 1)
#' y <- c(rnorm(200, 0, 1), rnorm(200, 8, 1))
#' init_binseg(y, X, gaussian_homo())
#' @export
init_binseg <- function(y, X, family, depth = 4L, min_seg = 10L) {
  n <- length(y)
  found <- integer(0)
  split_rec <- function(lo, hi, level) {
    if (level > depth || (hi - lo) < 2L * min_seg) return(invisible())
    idx <- (lo + 1L):hi
    res <- bruteforce(y[idx], X[idx, , drop = FALSE], family)
    b <- lo + res$bp
    if (b > lo + min_seg && b < hi - min_seg) {
      found <<- c(found, b)
      split_rec(lo, b, level + 1L)
      split_rec(b, hi, level + 1L)
    }
  }
  split_rec(0L, n, 1L)
  sort(unique(as.integer(found)))
}
