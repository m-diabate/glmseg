#' Test for a zero-length object
#' @param x an object
#' @return `TRUE` if `length(x) == 0`, else `FALSE`
#' @keywords internal
is_empty <- function(x) length(x) == 0L

#' Bidiagonal transition matrix (log scale)
#'
#' @param K number of segments (`K-1` breakpoints)
#' @return a `K x K` matrix in log scale: only "stay" (diagonal) and
#'   "advance by one" (super-diagonal) transitions are allowed.
#' @keywords internal
matPijk <- function(K) {
  pi_jk <- matrix(0, K, K)
  diag(pi_jk) <- 1
  if (K > 1L) {
    for (i in seq_len(K - 1L)) pi_jk[i, i + 1L] <- 1
  }
  log(pi_jk)
}

#' Thin a sorted vector so consecutive values differ by more than `tol`
#'
#' @param vect numeric vector
#' @param tol tolerance
#' @return a sorted subset of `vect`
#' @keywords internal
almost_unique_tol <- function(vect, tol) {
  if (length(vect) == 0L) return(vect)
  vect <- sort(vect)
  keep <- vect[1L]
  last <- vect[1L]
  for (v in vect[-1L]) {
    if (v - last > tol) {
      keep <- c(keep, v)
      last <- v
    }
  }
  keep
}
