#' max-Forward log-likelihood of the best segmentation
#'
#' Log-scaled max-forward recursion over a constrained HMM whose transition
#' matrix allows only "stay" and "advance by one".
#'
#' @param le n by K matrix of log-emission values.
#' @param lpi K by K log-transition matrix (see [matPijk()]).
#' @return scalar log-likelihood of the best segmentation.
#' @keywords internal
maxfw_loglik <- function(le, lpi) {
  d <- dim(le); n <- d[1L]; K <- d[2L]
  pi_jk <- exp(lpi)
  Fw <- matrix(NA_real_, n, K)
  L <- numeric(n)
  l <- apply(le, 1L, max)
  e <- exp(le - l)
  Fw[1L, ] <- 0; Fw[1L, 1L] <- e[1L, 1L]
  tmp <- max(Fw[1L, ]); L[1L] <- l[1L] + log(tmp); Fw[1L, ] <- Fw[1L, ] / tmp
  for (i in 2:n) {
    for (k in seq_len(K)) Fw[i, k] <- max(Fw[i - 1L, ] * pi_jk[, k] * e[i, k])
    tmp <- max(Fw[i, ]); L[i] <- L[i - 1L] + l[i] + log(tmp)
    Fw[i, ] <- Fw[i, ] / tmp
  }
  L[n] + log(Fw[n, K])
}

#' Best segmentation via max-Forward/max-Backward
#'
#' @inheritParams maxfw_loglik
#' @return list with `loglik`, `seg` (integer segment membership per point) and
#'   `bp` (breakpoint positions).
#' @keywords internal
maxfwbw_seg <- function(le, lpi) {
  d <- dim(le); n <- d[1L]; K <- d[2L]
  pi_jk <- exp(lpi)
  l <- apply(le, 1L, max)
  e <- exp(le - l)

  Fw <- matrix(NA_real_, n, K); L <- numeric(n)
  Fw[1L, ] <- 0; Fw[1L, 1L] <- e[1L, 1L]
  tmp <- max(Fw[1L, ]); L[1L] <- l[1L] + log(tmp); Fw[1L, ] <- Fw[1L, ] / tmp
  for (i in 2:n) {
    for (k in seq_len(K)) Fw[i, k] <- max(Fw[i - 1L, ] * pi_jk[, k] * e[i, k])
    tmp <- max(Fw[i, ]); L[i] <- L[i - 1L] + l[i] + log(tmp)
    Fw[i, ] <- Fw[i, ] / tmp
  }

  Bw <- matrix(NA_real_, n, K); M <- numeric(n)
  Bw[n, ] <- 0; Bw[n, K] <- 1
  tmp <- max(Bw[n, ]); M[n] <- log(tmp); Bw[n, ] <- Bw[n, ] / tmp
  for (i in n:2) {
    for (j in seq_len(K)) Bw[i - 1L, j] <- max(pi_jk[j, ] * e[i, ] * Bw[i, ])
    tmp <- max(Bw[i - 1L, ]); M[i - 1L] <- M[i] + l[i] + log(tmp)
    Bw[i - 1L, ] <- Bw[i - 1L, ] / tmp
  }

  post <- Fw * Bw
  seg <- apply(post, 1L, which.max)
  loglik <- L[n] + log(Fw[n, K])
  bp <- which(diff(seg) == 1L)
  names(bp) <- NULL
  list(loglik = loglik, seg = as.integer(seg), bp = as.integer(bp))
}
