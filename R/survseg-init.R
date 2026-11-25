# Binary-segmentation candidates for survival data.
survseg_init_binseg <- function(time, status, X, depth = 4L, min_seg = 15L) {
  n <- length(time)
  found <- integer(0)
  split_rec <- function(lo, hi, level) {
    if (level > depth || (hi - lo) < 2L * min_seg) return(invisible())
    idx <- (lo + 1L):hi
    res <- survseg_bruteforce(time[idx], status[idx], X[idx, , drop = FALSE])
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

# Fused-Lasso candidates for survival, via glmnet Cox on the triangular design.
survseg_init_fused_lasso <- function(time, status, K,
                                     min_seg = 50L, factor_lo = 5L) {
  n <- length(time)
  D <- matrix(1, n, n); D[upper.tri(D)] <- 0
  yy <- survival::Surv(time, status)
  # glmnet's Cox path can emit benign convergence warnings on the
  # over-parametrized triangular design; they do not affect candidate selection.
  fit <- suppressWarnings(glmnet::glmnet(x = D, y = yy, family = "cox"))
  n_seg <- apply(fit$beta, 2L, function(b) sum(abs(b) > 0))
  reached <- which(n_seg >= factor_lo * (K - 1L))
  if (length(reached) == 0L) {
    stop("Fused-Lasso (cox) initialization did not reach the requested ",
         "resolution; try init = 'binseg' or a different K.", call. = FALSE)
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
