#' Time-delay embedding
#'
#' Builds the Takens delay vectors \eqn{(x_i, x_{i+\tau}, \dots,
#' x_{i+(m-1)\tau})} of a scalar series.
#'
#' @param x Numeric vector.
#' @param m Embedding dimension.
#' @param tau Delay in samples.
#' @return A matrix with `length(x) - (m - 1) * tau` rows and `m` columns.
#' @examples
#' embed_delay(1:10, m = 3, tau = 2)
#' @export
embed_delay <- function(x, m = 3L, tau = 1L) {
  embed_cpp(as.numeric(x), as.integer(m), as.integer(tau))
}

# Turn whatever the user hands in into a matrix of state vectors.
as_states <- function(x, m, tau) {
  if (inherits(x, "chaos_trajectory")) {
    return(as.matrix(as.data.frame(unclass(x))[, c("x", "y", "z")]))
  }
  if (inherits(x, "chaos_map")) {
    x <- as.data.frame(unclass(x))
    return(as.matrix(x[, setdiff(names(x), "orbit")]))
  }
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.matrix(x)) {
    storage.mode(x) <- "double"
    return(x)
  }
  if (m > 1L) embed_delay(x, m, tau) else matrix(as.numeric(x), ncol = 1)
}

norm_code <- function(norm) {
  c(euclidean = 1L, max = 2L, manhattan = 3L)[[match.arg(norm, c("euclidean", "max", "manhattan"))]]
}

#' Choose a recurrence threshold for a target recurrence rate
#'
#' Returns the distance quantile that makes roughly a fraction `rr` of all
#' state pairs recurrent. For large `N` the pairwise distances are sampled.
#'
#' @param X Matrix of state vectors (rows).
#' @param rr Target recurrence rate in `(0, 1)`.
#' @param norm `"euclidean"`, `"max"` or `"manhattan"`.
#' @param n_pairs Maximum number of sampled pairs.
#' @return A single threshold value.
#' @export
choose_eps <- function(X, rr = 0.05, norm = "euclidean", n_pairs = 1e6) {
  d <- sample_distances_cpp(X, as.integer(n_pairs), norm_code(norm))
  unname(quantile(d, rr, names = FALSE))
}

#' Recurrence matrix
#'
#' Computes \eqn{R_{ij} = 1} if \eqn{\|x_i - x_j\| \le \epsilon}. The input can
#' be a scalar series (delay-embedded with `m` and `tau`), a matrix of state
#' vectors, or a `chaos_trajectory` / `chaos_map` object.
#'
#' @param x Series, matrix or chaosrcpp object.
#' @param eps Threshold. When `NULL` it is chosen with [choose_eps()] so that
#'   the recurrence rate is about `rr`.
#' @param rr Target recurrence rate used when `eps` is `NULL`.
#' @param m,tau Embedding dimension and delay for scalar input.
#' @param norm Distance norm.
#' @return A logical matrix with attributes `eps`, `norm`, and class
#'   `recurrence_matrix`.
#' @examples
#' x <- sin(seq(0, 20 * pi, length.out = 600))
#' R <- recurrence_matrix(x, m = 2, tau = 8, rr = 0.1)
#' plot_recurrence(R)
#' @seealso [rqa()], [plot_recurrence()]
#' @export
recurrence_matrix <- function(x, eps = NULL, rr = 0.05, m = 1L, tau = 1L,
                              norm = c("euclidean", "max", "manhattan")) {
  norm <- match.arg(norm)
  X <- as_states(x, m, tau)
  if (is.null(eps)) eps <- choose_eps(X, rr, norm)
  R <- recurrence_matrix_cpp(X, eps, norm_code(norm))
  attr(R, "eps") <- eps
  attr(R, "norm") <- norm
  class(R) <- c("recurrence_matrix", class(R))
  R
}

#' Recurrence Quantification Analysis
#'
#' Computes the standard RQA measures of a recurrence matrix in one pass
#' through compiled code: recurrence rate (`RR`), determinism (`DET`),
#' average and maximal diagonal line length (`L`, `Lmax`), divergence
#' (`DIV = 1/Lmax`), Shannon entropy of the diagonal line length distribution
#' (`ENTR`), laminarity (`LAM`), trapping time (`TT`), longest vertical line
#' (`Vmax`) and `RATIO = DET / RR`. Definitions follow Marwan et al. (2007).
#'
#' Memory use is `N^2` bytes for `N` state vectors, so `N` up to about ten
#' thousand is comfortable; thin or window longer series.
#'
#' @inheritParams recurrence_matrix
#' @param lmin Minimal diagonal line length counted in `DET`, `L`, `ENTR`.
#' @param vmin Minimal vertical line length counted in `LAM`, `TT`.
#' @param theiler Theiler window: cells with `|i - j| < theiler` are ignored.
#'   `1` excludes only the main diagonal.
#' @return A list of class `rqa` with the measures, the settings and the
#'   diagonal / vertical line length histograms.
#' @references Marwan, N., Romano, M. C., Thiel, M., & Kurths, J. (2007).
#'   Recurrence plots for the analysis of complex systems. Physics Reports,
#'   438(5-6), 237-329.
#' @examples
#' set.seed(1)
#' periodic <- sin(seq(0, 30 * pi, length.out = 1500))
#' noise <- rnorm(1500)
#' rqa(periodic, m = 3, tau = 10)
#' rqa(noise, m = 3, tau = 1)
#' @seealso [rqa_features()], [recurrence_matrix()]
#' @export
rqa <- function(x, eps = NULL, rr = 0.05, m = 1L, tau = 1L,
                norm = c("euclidean", "max", "manhattan"),
                lmin = 2L, vmin = 2L, theiler = 1L) {
  norm <- match.arg(norm)
  X <- as_states(x, m, tau)
  if (nrow(X) > 20000L) {
    warning(sprintf("%d state vectors need %.1f GB for the recurrence matrix; consider thinning.",
                    nrow(X), nrow(X)^2 / 1e9))
  }
  if (is.null(eps)) eps <- choose_eps(X, rr, norm)
  res <- rqa_cpp(X, eps, norm_code(norm), as.integer(lmin), as.integer(vmin),
                 as.integer(theiler))
  res$settings <- list(eps = eps, norm = norm, m = m, tau = tau, lmin = lmin,
                       vmin = vmin, theiler = theiler)
  class(res) <- "rqa"
  res
}

rqa_measure_names <- c("RR", "DET", "L", "Lmax", "DIV", "ENTR", "LAM", "TT", "Vmax", "RATIO")

#' RQA measures as a one-row data frame
#'
#' Convenience wrapper around [rqa()] that returns only the ten measures, in
#' the shape needed for feature tables (one row per signal or window).
#'
#' @inheritParams rqa
#' @param ... Passed to [rqa()].
#' @return A one-row data frame.
#' @examples
#' rqa_features(sin(seq(0, 30 * pi, length.out = 1000)), m = 2, tau = 5)
#' @export
rqa_features <- function(x, ...) {
  r <- rqa(x, ...)
  as.data.frame(r[rqa_measure_names])
}

#' @export
as.data.frame.rqa <- function(x, ...) as.data.frame(unclass(x)[rqa_measure_names])

#' @export
print.rqa <- function(x, digits = 4, ...) {
  s <- x$settings
  cat(sprintf("Recurrence Quantification Analysis (%d states, eps = %s, %s norm, m = %d, tau = %d, theiler = %d)\n",
              x$n_states, format(s$eps, digits = 4), s$norm, s$m, s$tau, s$theiler))
  vals <- unlist(unclass(x)[rqa_measure_names])
  print(round(vals, digits))
  invisible(x)
}

#' Draw a recurrence plot
#'
#' Base-graphics raster of a recurrence matrix with the conventional
#' orientation (time along both axes, origin bottom-left).
#'
#' @param R A logical matrix from [recurrence_matrix()], or anything that
#'   [recurrence_matrix()] accepts (then `...` is passed on to it).
#' @param col Colour of recurrence points.
#' @param main Title.
#' @param ... Extra arguments for [recurrence_matrix()] when `R` is not
#'   already a matrix.
#' @return Invisibly, the recurrence matrix.
#' @examples
#' tr <- lorenz(n = 3000, dt = 0.02)
#' plot_recurrence(tr, rr = 0.05)
#' @export
plot_recurrence <- function(R, col = "#1f2a44", main = "Recurrence plot", ...) {
  if (!is.matrix(R)) R <- recurrence_matrix(R, ...)
  n <- nrow(R)
  op <- par(mar = c(4, 4, 2.5, 1), pty = "s")
  on.exit(par(op))
  image(seq_len(n), seq_len(n), unclass(R) * 1, col = c("white", col), useRaster = TRUE,
        xlab = "i", ylab = "j", main = main, axes = FALSE)
  axis(1); axis(2); box()
  invisible(R)
}
