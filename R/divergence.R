#' The butterfly effect of rounding error
#'
#' Integrates the Lorenz system from one initial condition in several
#' floating-point precisions and measures how far each trajectory drifts from
#' a high-precision reference. Since the algorithm, the step size and the
#' start point are identical, the only source of disagreement is rounding.
#' On the chaotic attractor that disagreement grows like
#' \eqn{\epsilon \, e^{\lambda_1 t}}, where \eqn{\epsilon} is the unit
#' roundoff (about 6e-8 for fp32, 1e-16 for fp64) and \eqn{\lambda_1 \approx
#' 0.906} is the largest Lyapunov exponent. Single precision therefore loses
#' the trajectory roughly twice as early as double precision, and both lose it
#' eventually.
#'
#' The default reference is the compiled double-double level (`"dd"`, unit
#' roundoff about 5e-32), which is fast. `reference = "mpfr"` uses Rmpfr as the
#' gold standard instead; it is an R-level loop running at roughly a hundred
#' steps per second, so budget a minute for the default `t_max = 60`. With an
#' mpfr reference the dd trajectory is included in the comparison, so you get
#' three curves diverging at roughly `t = 20`, `40` and `80`.
#'
#' @inheritParams lorenz
#' @param t_max Length of the integration in time units.
#' @param precisions Precisions to compare against the reference. Defaults
#'   to every compiled level other than the reference.
#' @param reference `"dd"` (default), `"mpfr"` (needs Rmpfr) or `"fp64"`.
#' @param mpfr_bits Mantissa bits of the mpfr reference.
#' @return An object of class `precision_divergence`: a list with
#'   `trajectories` (named list of `chaos_trajectory` data frames including
#'   the reference), `distance` (long data frame with `t`, `precision`,
#'   `distance`), `divergence_time` (first time the distance exceeds
#'   `threshold`), `reference` and `params`.
#' @param threshold Distance at which a trajectory is declared lost.
#' @examples
#' \donttest{
#' d <- precision_divergence(t_max = 50)
#' d
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot_divergence(d)
#' }
#' @seealso [plot_divergence()]
#' @export
precision_divergence <- function(x0 = c(1, 1, 1), sigma = 10, rho = 28,
                                 beta = 8 / 3, dt = 0.01, t_max = 60,
                                 precisions = c("fp32", "fp64", "dd"),
                                 reference = c("dd", "mpfr", "fp64"),
                                 mpfr_bits = 128L, thin = 1L, threshold = 1) {
  reference <- match.arg(reference)
  precisions <- match.arg(precisions, c("fp32", "fp64", "dd"), several.ok = TRUE)
  precisions <- setdiff(precisions, reference)
  if (!length(precisions)) stop("nothing left to compare with the reference")
  n <- as.integer(round(t_max / dt))

  run <- function(prec) {
    lorenz(x0 = x0, sigma = sigma, rho = rho, beta = beta, dt = dt, n = n,
           precision = prec, thin = thin, mpfr_bits = mpfr_bits)
  }
  ref <- run(reference)
  trajs <- lapply(precisions, run)
  names(trajs) <- precisions

  dist <- do.call(rbind, lapply(precisions, function(prec) {
    data.frame(t = trajs[[prec]]$t, precision = prec,
               distance = trajectory_distance(trajs[[prec]], ref))
  }))
  dist$precision <- factor(dist$precision, levels = precisions)

  div_time <- vapply(precisions, function(prec) {
    dd <- dist[dist$precision == prec, ]
    idx <- which(dd$distance > threshold)
    if (length(idx)) dd$t[idx[1]] else NA_real_
  }, numeric(1))

  trajs[[reference]] <- ref
  structure(list(trajectories = trajs, distance = dist,
                 divergence_time = div_time, reference = reference,
                 threshold = threshold, dt = dt, t_max = t_max,
                 params = list(sigma = sigma, rho = rho, beta = beta, x0 = x0)),
            class = "precision_divergence")
}

# Euclidean distance between two trajectories, computed in high precision
# when either side carries more digits than a double (dd or mpfr).
trajectory_distance <- function(a, b) {
  high <- c("dd", "mpfr")
  if ((attr(a, "precision") %in% high || attr(b, "precision") %in% high) &&
      requireNamespace("Rmpfr", quietly = TRUE)) {
    bits <- max(128L, if (!is.null(attr(a, "mpfr"))) Rmpfr::getPrec(attr(a, "mpfr"))[1] else 0L,
                if (!is.null(attr(b, "mpfr"))) Rmpfr::getPrec(attr(b, "mpfr"))[1] else 0L)
    d <- exact_states(a, bits) - exact_states(b, bits)
    return(as.numeric(sqrt(d[, 1]^2 + d[, 2]^2 + d[, 3]^2)))
  }
  sqrt((a$x - b$x)^2 + (a$y - b$y)^2 + (a$z - b$z)^2)
}

#' @export
print.precision_divergence <- function(x, ...) {
  cat(sprintf("<precision_divergence> Lorenz, t in [0, %g], dt = %g, reference = %s\n",
              x$t_max, x$dt, x$reference))
  cat(sprintf("Time until distance to reference exceeds %g:\n", x$threshold))
  for (p in names(x$divergence_time)) {
    v <- x$divergence_time[[p]]
    cat(sprintf("  %-5s %s\n", p, if (is.na(v)) "never (within t_max)" else sprintf("t = %.2f", v)))
  }
  invisible(x)
}

#' Plot a precision-divergence experiment
#'
#' Two stacked panels: the `x(t)` traces of every precision including the
#' reference, and the distance to the reference on a log scale. Dashed lines
#' show the textbook prediction \eqn{\epsilon \, e^{\lambda_1 t}} for each
#' precision.
#'
#' @param x A `precision_divergence` object.
#' @param lambda1 Largest Lyapunov exponent used for the dashed prediction
#'   lines (0.9056 for the classic Lorenz parameters). Set to `NULL` to hide.
#' @param trace_var Which coordinate to draw in the top panel.
#' @param ... Unused.
#' @return A ggplot / patchwork object (patchwork optional; without it a list
#'   of two ggplots is returned).
#' @export
plot_divergence <- function(x, lambda1 = 0.9056, trace_var = c("x", "y", "z"), ...) {
  need_pkg("ggplot2", "plot_divergence()")
  trace_var <- match.arg(trace_var)
  traces <- do.call(rbind, lapply(names(x$trajectories), function(p) {
    tr <- x$trajectories[[p]]
    data.frame(t = tr$t, value = tr[[trace_var]], precision = p)
  }))
  lv <- c(setdiff(names(x$trajectories), x$reference), x$reference)
  traces$precision <- factor(traces$precision, levels = lv)
  cols <- c(fp32 = "#D55E00", fp64 = "#0072B2", dd = "#009E73", mpfr = "#333333")

  p1 <- ggplot2::ggplot(traces, ggplot2::aes(t, value, colour = precision)) +
    ggplot2::geom_line(linewidth = 0.4, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::labs(y = paste0(trace_var, "(t)"), x = NULL, colour = "precision",
                  title = "Same initial condition, same algorithm, same step size",
                  subtitle = "Only the floating-point precision differs") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top")

  dist <- x$distance
  dist <- dist[dist$distance > 0, ]
  p2 <- ggplot2::ggplot(dist, ggplot2::aes(t, distance, colour = precision)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = x$threshold, linetype = "dotted", colour = "grey40") +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    ggplot2::scale_y_log10(labels = function(v) sprintf("1e%+d", round(log10(v)))) +
    ggplot2::labs(y = sprintf("distance to %s reference", x$reference), x = "t",
                  caption = if (!is.null(lambda1)) "dashed: unit roundoff times exp(lambda1 t)" else NULL) +
    ggplot2::theme_minimal(base_size = 12)

  if (!is.null(lambda1)) {
    tt <- seq(0, x$t_max, length.out = 200)
    theory <- do.call(rbind, lapply(levels(dist$precision), function(p) {
      eps <- precision_epsilon[[p]]
      data.frame(t = tt, distance = pmin(eps * exp(lambda1 * tt), 100), precision = p)
    }))
    theory$precision <- factor(theory$precision, levels = levels(dist$precision))
    p2 <- p2 + ggplot2::geom_line(data = theory, linetype = "dashed", linewidth = 0.5, alpha = 0.7)
  }

  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 1.2))
  } else {
    list(traces = p1, distance = p2)
  }
}
