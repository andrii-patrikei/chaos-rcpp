#' Logistic map orbit
#'
#' Iterates \eqn{x_{n+1} = r x_n (1 - x_n)}.
#'
#' @param r Growth parameter in `[0, 4]`.
#' @param x0 Initial value in `(0, 1)`.
#' @param n Number of iterations.
#' @return Numeric vector `x_0, ..., x_n`.
#' @examples
#' plot(logistic_orbit(3.9, n = 100), type = "b", pch = 16, cex = 0.6)
#' @export
logistic_orbit <- function(r, x0 = 0.2, n = 100L) {
  logistic_orbit_cpp(x0, r, as.integer(n))
}

#' Bifurcation diagram of the logistic map
#'
#' For every value of `r` the map is iterated `n_transient` times to forget
#' the initial condition and the next `n_keep` iterates are recorded. With
#' the defaults this is 400,000 recorded points from 1.4 million iterations,
#' which the compiled loop does in a fraction of a second.
#'
#' @param r Vector of parameter values.
#' @param x0 Initial value. Avoid `0.5`: it is the critical point of the map,
#'   so at \eqn{r = 4} it lands on the fixed point 0 after two steps.
#' @param n_transient Iterates discarded per `r`.
#' @param n_keep Iterates recorded per `r`.
#' @return Data frame with columns `r` and `x`, class `bifurcation`.
#' @examples
#' bif <- logistic_bifurcation(seq(2.8, 4, length.out = 800))
#' plot(bif)
#' @seealso [logistic_lyapunov()], [plot_bifurcation()]
#' @export
logistic_bifurcation <- function(r = seq(2.5, 4, length.out = 2000), x0 = 0.1,
                                 n_transient = 500L, n_keep = 200L) {
  out <- logistic_bifurcation_cpp(as.numeric(r), x0, as.integer(n_transient),
                                  as.integer(n_keep))
  df <- data.frame(r = out$r, x = out$x)
  attr(df, "map") <- "logistic"
  attr(df, "parameter") <- "r"
  class(df) <- c("bifurcation", "data.frame")
  df
}

#' Lyapunov exponent of the logistic map
#'
#' \eqn{\lambda(r) = \lim_{n\to\infty} \frac{1}{n} \sum_k \log |r (1 - 2 x_k)|}.
#' The exponent is negative on periodic windows, zero at bifurcation points
#' and positive in the chaotic regime; at \eqn{r = 4} it equals \eqn{\log 2}.
#'
#' @inheritParams logistic_bifurcation
#' @param n Iterates averaged per `r`.
#' @return Data frame with columns `r` and `lambda`.
#' @examples
#' ly <- logistic_lyapunov(c(3.2, 3.5, 3.9, 4))
#' ly
#' @export
logistic_lyapunov <- function(r = seq(2.5, 4, length.out = 2000), x0 = 0.1,
                              n_transient = 500L, n = 2000L) {
  lambda <- logistic_lyapunov_cpp(as.numeric(r), x0, as.integer(n_transient),
                                  as.integer(n))
  data.frame(r = r, lambda = lambda)
}

#' Cobweb (return-map) construction for the logistic map
#'
#' Produces the line segments of the classic graphical iteration: from
#' `(x_n, x_n)` up to the curve `(x_n, f(x_n))` and across to the diagonal
#' `(f(x_n), f(x_n))`.
#'
#' @inheritParams logistic_orbit
#' @return A list with the orbit and a data frame of segments
#'   (`x`, `y`, `xend`, `yend`, `step`), class `cobweb`.
#' @examples
#' cw <- logistic_cobweb(3.7, x0 = 0.1, n = 60)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot_cobweb(cw)
#' @export
logistic_cobweb <- function(r, x0 = 0.1, n = 50L) {
  x <- logistic_orbit_cpp(x0, r, as.integer(n))
  k <- seq_len(n)
  up <- data.frame(x = x[k], y = x[k], xend = x[k], yend = x[k + 1], step = k)
  across <- data.frame(x = x[k], y = x[k + 1], xend = x[k + 1], yend = x[k + 1], step = k)
  seg <- rbind(up, across)
  seg <- seg[order(seg$step), ]
  # the very first segment starts on the x axis
  seg$y[1] <- 0
  structure(list(r = r, orbit = x, segments = seg), class = "cobweb")
}

#' Henon map attractor
#'
#' Iterates \eqn{x_{n+1} = 1 - a x_n^2 + y_n, \; y_{n+1} = b x_n}. The classic
#' strange attractor appears at `a = 1.4, b = 0.3`.
#'
#' @param a,b Map parameters.
#' @param x0,y0 Initial condition.
#' @param n Points returned.
#' @param n_transient Iterates discarded first.
#' @return Data frame with columns `x`, `y`, class `chaos_map`.
#' @examples
#' plot(henon(n = 20000))
#' @export
henon <- function(a = 1.4, b = 0.3, x0 = 0, y0 = 0, n = 100000L, n_transient = 100L) {
  m <- henon_cpp(x0, y0, a, b, as.integer(n), as.integer(n_transient))
  new_map(m, "henon", list(a = a, b = b))
}

#' Bifurcation diagram of the Henon map (in `a`, for fixed `b`)
#'
#' @inheritParams henon
#' @param a Vector of `a` values.
#' @param b Fixed value of the second parameter.
#' @param n_keep Iterates recorded per `a`.
#' @return Data frame with columns `a` and `x`, class `bifurcation`.
#' @examples
#' plot(henon_bifurcation(seq(1, 1.42, length.out = 600)))
#' @export
henon_bifurcation <- function(a = seq(1, 1.42, length.out = 1500), b = 0.3,
                              x0 = 0, y0 = 0, n_transient = 500L, n_keep = 150L) {
  out <- henon_bifurcation_cpp(as.numeric(a), b, x0, y0, as.integer(n_transient),
                               as.integer(n_keep))
  df <- data.frame(a = out$a, x = out$x)
  attr(df, "map") <- "henon"
  attr(df, "parameter") <- "a"
  class(df) <- c("bifurcation", "data.frame")
  df
}

#' Ikeda map attractor
#'
#' A model of light in a nonlinear optical cavity:
#' \eqn{t_n = 0.4 - 6 / (1 + x_n^2 + y_n^2)},
#' \eqn{x_{n+1} = 1 + u (x_n \cos t_n - y_n \sin t_n)},
#' \eqn{y_{n+1} = u (x_n \sin t_n + y_n \cos t_n)}. The attractor for
#' `u = 0.918` is one of the prettiest in the business.
#'
#' @param u Dissipation parameter, chaotic for `u >= 0.6` or so.
#' @inheritParams henon
#' @return Data frame with columns `x`, `y`, class `chaos_map`.
#' @examples
#' plot(ikeda(n = 50000))
#' @export
ikeda <- function(u = 0.918, x0 = 0.1, y0 = 0.1, n = 200000L, n_transient = 100L) {
  m <- ikeda_cpp(x0, y0, u, as.integer(n), as.integer(n_transient))
  new_map(m, "ikeda", list(u = u))
}

#' Chirikov standard map
#'
#' Iterates \eqn{p_{n+1} = p_n + K \sin \theta_n}, \eqn{\theta_{n+1} = \theta_n +
#' p_{n+1}} (both modulo \eqn{2\pi}) for many initial conditions at once, so
#' that the phase portrait shows KAM tori, resonance islands and the chaotic
#' sea side by side. `K = 0.971635` is the critical value at which the last
#' invariant torus breaks and global chaos sets in.
#'
#' @param K Kick strength.
#' @param n_orbits Number of random initial conditions.
#' @param n Iterates per orbit.
#' @param seed Seed for the initial conditions (`NULL` to leave the RNG alone).
#' @return Data frame with columns `orbit`, `theta`, `p`, class `chaos_map`.
#' @examples
#' plot(standard_map(K = 0.6, n_orbits = 40, n = 400))
#' @export
standard_map <- function(K = 0.971635, n_orbits = 60L, n = 1000L, seed = 1L) {
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
  }
  theta0 <- stats::runif(n_orbits, 0, 2 * pi)
  p0 <- stats::runif(n_orbits, 0, 2 * pi)
  m <- standard_map_cpp(theta0, p0, K, as.integer(n))
  df <- data.frame(orbit = as.integer(m[, 1]), theta = m[, 2], p = m[, 3])
  attr(df, "map") <- "standard"
  attr(df, "params") <- list(K = K)
  class(df) <- c("chaos_map", "data.frame")
  df
}

new_map <- function(m, map, params) {
  df <- data.frame(x = m[, 1], y = m[, 2])
  attr(df, "map") <- map
  attr(df, "params") <- params
  class(df) <- c("chaos_map", "data.frame")
  df
}

#' @export
print.chaos_map <- function(x, ...) {
  p <- attr(x, "params")
  cat(sprintf("<chaos_map> %s map, %d points\n", attr(x, "map"), nrow(x)))
  cat("  parameters:", paste(names(p), signif(unlist(p), 6), sep = " = ", collapse = ", "), "\n")
  print(utils::head(as.data.frame(unclass(x)), 6), ...)
  invisible(x)
}

# --------------------------------------------------------------------------
# Base-graphics plot methods (fast, dependency-free)
# --------------------------------------------------------------------------

#' Plot a map attractor or a bifurcation diagram
#'
#' Base-graphics methods that draw hundreds of thousands of points quickly.
#' For publication-style output see [plot_bifurcation()] and
#' [plot_standard_map()] (ggplot2).
#'
#' @param x A `chaos_map` or `bifurcation` object.
#' @param col Point colour.
#' @param pch,cex Point symbol and size; the default is a single pixel.
#' @param ... Passed on to [graphics::plot()].
#' @return Invisibly, `x`.
#' @export
plot.chaos_map <- function(x, col = "#1f2a44", pch = ".", cex = 1, ...) {
  map <- attr(x, "map")
  op <- par(mar = c(4, 4, 2.5, 1))
  on.exit(par(op))
  if (identical(map, "standard")) {
    pal <- hcl.colors(max(x$orbit), "Zissou 1")
    plot(x$theta, x$p, col = pal[x$orbit], pch = pch, cex = cex, asp = 1,
         xlab = expression(theta), ylab = "p", xaxs = "i", yaxs = "i",
         main = sprintf("Standard map, K = %g", attr(x, "params")$K), ...)
  } else {
    plot(x$x, x$y, col = col, pch = pch, cex = cex, xlab = "x", ylab = "y",
         main = paste(tools::toTitleCase(map), "map"), ...)
  }
  invisible(x)
}

#' @rdname plot.chaos_map
#' @export
plot.bifurcation <- function(x, col = "#1f2a44", pch = ".", cex = 1, ...) {
  prm <- attr(x, "parameter")
  op <- par(mar = c(4, 4, 2.5, 1))
  on.exit(par(op))
  plot(x[[prm]], x$x, col = adjustcolor(col, 0.5), pch = pch, cex = cex,
       xlab = prm, ylab = "x", xaxs = "i",
       main = paste(tools::toTitleCase(attr(x, "map")), "map bifurcation diagram"), ...)
  invisible(x)
}

# --------------------------------------------------------------------------
# ggplot2 helpers
# --------------------------------------------------------------------------

#' Bifurcation diagram with the Lyapunov exponent underneath
#'
#' The upper panel is the bifurcation diagram, the lower one the Lyapunov
#' exponent over the same parameter range with the chaotic region
#' (\eqn{\lambda > 0}) shaded. Lining them up shows that chaos begins exactly
#' where the exponent crosses zero and that every periodic window is a dip
#' below zero.
#'
#' @param bif A `bifurcation` data frame, see [logistic_bifurcation()].
#' @param lyap Optional data frame from [logistic_lyapunov()] over the same
#'   parameter range. When `NULL` only the diagram is drawn.
#' @param point_alpha Transparency of the bifurcation points.
#' @param colour Point colour.
#' @return A patchwork object when `lyap` is given and patchwork is installed,
#'   otherwise a ggplot (or a list of two ggplots).
#' @examples
#' \donttest{
#' r <- seq(2.8, 4, length.out = 1500)
#' plot_bifurcation(logistic_bifurcation(r), logistic_lyapunov(r))
#' }
#' @export
plot_bifurcation <- function(bif, lyap = NULL, point_alpha = 0.08, colour = "#1f2a44") {
  need_pkg("ggplot2", "plot_bifurcation()")
  prm <- attr(bif, "parameter")
  dat <- data.frame(p = bif[[prm]], x = bif$x)
  p1 <- ggplot2::ggplot(dat, ggplot2::aes(p, x)) +
    ggplot2::geom_point(shape = ".", alpha = point_alpha, colour = colour) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = if (is.null(lyap)) prm else NULL, y = "x",
                  title = paste(tools::toTitleCase(attr(bif, "map")), "map")) +
    ggplot2::theme_minimal(base_size = 12)
  if (is.null(lyap)) return(p1)

  ly <- data.frame(p = lyap[[1]], lambda = lyap$lambda)
  ly$sign <- ifelse(ly$lambda > 0, "chaotic", "regular")
  p2 <- ggplot2::ggplot(ly, ggplot2::aes(p, lambda)) +
    ggplot2::geom_area(data = transform(ly, lambda = pmax(lambda, 0)),
                       fill = "#D55E00", alpha = 0.25) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    ggplot2::geom_line(colour = colour, linewidth = 0.35) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(max(-4, min(ly$lambda, na.rm = TRUE)), max(ly$lambda, na.rm = TRUE) * 1.05)) +
    ggplot2::labs(x = prm, y = expression(lambda), caption = "shaded: lambda > 0 (chaos)") +
    ggplot2::theme_minimal(base_size = 12)
  p1 <- p1 + ggplot2::theme(axis.text.x = ggplot2::element_blank())

  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(2, 1))
  } else {
    list(bifurcation = p1, lyapunov = p2)
  }
}

#' Plot a cobweb diagram
#'
#' @param cw A `cobweb` object from [logistic_cobweb()].
#' @param ... Unused.
#' @return A ggplot.
#' @export
plot_cobweb <- function(cw, ...) {
  need_pkg("ggplot2", "plot_cobweb()")
  xs <- seq(0, 1, length.out = 400)
  curve <- data.frame(x = xs, y = cw$r * xs * (1 - xs))
  seg <- cw$segments
  ggplot2::ggplot() +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_line(data = curve, ggplot2::aes(x, y), colour = "#1f2a44", linewidth = 0.7) +
    ggplot2::geom_segment(data = seg,
                          ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = step),
                          linewidth = 0.4) +
    ggplot2::scale_colour_viridis_c(option = "C", end = 0.9, guide = "none") +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = expression(x[n]), y = expression(x[n + 1]),
                  title = sprintf("Logistic map cobweb, r = %g, x0 = %g", cw$r, cw$orbit[1])) +
    ggplot2::theme_minimal(base_size = 12)
}

#' Plot the standard-map phase portrait with ggplot2
#'
#' @param sm Output of [standard_map()].
#' @param size Point size.
#' @return A ggplot.
#' @export
plot_standard_map <- function(sm, size = 0.15) {
  need_pkg("ggplot2", "plot_standard_map()")
  ggplot2::ggplot(sm, ggplot2::aes(theta, p, colour = factor(orbit))) +
    ggplot2::geom_point(size = size, alpha = 0.8) +
    ggplot2::scale_colour_manual(values = hcl.colors(max(sm$orbit), "Zissou 1"), guide = "none") +
    ggplot2::coord_equal(xlim = c(0, 2 * pi), ylim = c(0, 2 * pi), expand = FALSE) +
    ggplot2::labs(x = expression(theta), y = "p",
                  title = sprintf("Chirikov standard map, K = %g", attr(sm, "params")$K),
                  subtitle = "KAM tori, resonance islands and the chaotic sea") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#0b0e1a", colour = NA),
                   panel.grid = ggplot2::element_blank())
}
