# Simulate right-censored Weibull-AFT data with K segments (intercept-only or
# with covariates). par = c(betas [p*K], shapes [K]); cens_scale per segment.
sim_survseg <- function(X, par, bp = NULL, cens_scale) {
  n <- nrow(X); p <- ncol(X); K <- length(bp) + 1L
  beta <- matrix(par[seq_len(p * K)], nrow = p)
  Xbeta <- X %*% beta
  bp0 <- c(0L, bp, n)
  time <- numeric(n); status <- numeric(n)
  for (k in seq_len(K)) {
    idx <- (bp0[k] + 1L):bp0[k + 1L]
    shape_k <- par[p * K + k]
    scale_k <- exp(Xbeta[idx, k])
    tt <- stats::rweibull(length(idx), shape = shape_k, scale = scale_k)
    cc <- stats::rexp(length(idx), rate = 1 / cens_scale[k])
    time[idx] <- pmin(tt, cc)
    status[idx] <- as.numeric(tt <= cc)
  }
  list(time = time, status = status)
}
