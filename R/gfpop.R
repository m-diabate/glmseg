#' Build a K-1 breakpoint chain graph for gfpop
#'
#' @param n_bp number of breakpoints (K-1).
#' @return a `gfpop` graph object.
#' @keywords internal
my_graph_fun <- function(n_bp) {
  edges <- lapply(seq_len(n_bp), function(i)
    gfpop::Edge(as.character(i), as.character(i + 1L), type = "std", penalty = 0))
  args <- c(edges,
            list(gfpop::StartEnd(start = "1", end = as.character(n_bp + 1L)),
                 all.null.edges = TRUE))
  do.call(gfpop::graph, args)
}

#' GFPOP mean-model segmentation (comparison method)
#'
#' @param y response vector.
#' @param K number of segments.
#' @return list with `bp` (breakpoints) and the raw gfpop `fit`.
#' @examples
#' if (requireNamespace("gfpop", quietly = TRUE)) {
#'   set.seed(1)
#'   y <- c(rnorm(200, 0, 1), rnorm(200, 8, 1))
#'   gfpop_gauss(y, K = 2)$bp
#' }
#' @export
gfpop_gauss <- function(y, K) {
  if (!requireNamespace("gfpop", quietly = TRUE)) {
    stop("Package 'gfpop' is required for gfpop_gauss(); please install it.",
         call. = FALSE)
  }
  g <- my_graph_fun(K - 1L)
  fit <- gfpop::gfpop(data = y, mygraph = g, type = "mean")
  list(bp = fit$changepoints[seq_len(K - 1L)], fit = fit)
}
