#' Plot a trajectory as a 2D phase portrait
#'
#' @param x A `chaos_trajectory` (from [lorenz()] or [roessler()]).
#' @param vars Two variable names to draw, e.g. `c("x", "z")`.
#' @param col Line colour, or a vector of colours (one per row) to encode time.
#' @param lwd Line width.
#' @param ... Passed on to [graphics::plot()].
#' @return Invisibly, `x`.
#' @examples
#' plot(lorenz(n = 8000))
#' @export
plot_attractor <- function(x, vars = c("x", "z"), col = NULL, lwd = 0.5, ...) {
  stopifnot(length(vars) == 2, all(vars %in% c("x", "y", "z")))
  op <- par(mar = c(4, 4, 2.5, 1))
  on.exit(par(op))
  ttl <- sprintf("%s attractor (%s)", tools::toTitleCase(attr(x, "system")), attr(x, "precision"))
  if (is.null(col)) {
    plot(x[[vars[1]]], x[[vars[2]]], type = "l", col = "#1f2a44", lwd = lwd,
         xlab = vars[1], ylab = vars[2], main = ttl, ...)
  } else {
    n <- nrow(x)
    plot(x[[vars[1]]], x[[vars[2]]], type = "n", xlab = vars[1], ylab = vars[2], main = ttl, ...)
    if (length(col) == 1L) col <- rep(col, n)
    graphics::segments(x[[vars[1]]][-n], x[[vars[2]]][-n], x[[vars[1]]][-1], x[[vars[2]]][-1],
                       col = col[-n], lwd = lwd)
  }
  invisible(x)
}

#' @rdname plot_attractor
#' @export
plot.chaos_trajectory <- function(x, vars = c("x", "z"), col = NULL, lwd = 0.5, ...) {
  plot_attractor(x, vars = vars, col = col, lwd = lwd, ...)
}

#' Interactive 3D attractor (plotly)
#'
#' Rotatable, zoomable 3D line plot of a trajectory, coloured by time or by
#' any column. Several trajectories can be overlaid to see, for instance, the
#' fp32 and fp64 runs of the same initial condition part ways.
#'
#' @param ... One or more `chaos_trajectory` objects. Names are used in the
#'   legend; unnamed ones are labelled by their precision.
#' @param color_by Column used for colour when a single trajectory is given
#'   (`"t"` by default).
#' @param width Line width.
#' @return A plotly htmlwidget.
#' @examples
#' \donttest{
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   plot_attractor_3d(lorenz(n = 6000))
#'   plot_attractor_3d(fp32 = lorenz(n = 3000, precision = "fp32"),
#'                     fp64 = lorenz(n = 3000))
#' }
#' }
#' @export
plot_attractor_3d <- function(..., color_by = "t", width = 2) {
  need_pkg("plotly", "plot_attractor_3d()")
  trajs <- list(...)
  if (length(trajs) == 1L && !inherits(trajs[[1]], "chaos_trajectory") && is.list(trajs[[1]])) {
    trajs <- trajs[[1]]
  }
  nm <- names(trajs)
  if (is.null(nm)) nm <- rep("", length(trajs))
  for (i in seq_along(trajs)) if (nm[i] == "") nm[i] <- attr(trajs[[i]], "precision") %||% paste0("run ", i)

  p <- plotly::plot_ly()
  if (length(trajs) == 1L) {
    tr <- trajs[[1]]
    p <- plotly::add_trace(p, data = as.data.frame(unclass(tr)), x = ~x, y = ~y, z = ~z,
                           type = "scatter3d", mode = "lines",
                           line = list(width = width, color = tr[[color_by]],
                                       colorscale = "Viridis", colorbar = list(title = color_by)),
                           name = nm[1], hoverinfo = "x+y+z")
  } else {
    cols <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#333333")
    for (i in seq_along(trajs)) {
      tr <- trajs[[i]]
      p <- plotly::add_trace(p, data = as.data.frame(unclass(tr)), x = ~x, y = ~y, z = ~z,
                             type = "scatter3d", mode = "lines",
                             line = list(width = width, color = cols[(i - 1) %% length(cols) + 1]),
                             name = nm[i], hoverinfo = "x+y+z+name")
    }
  }
  sysname <- attr(trajs[[1]], "system") %||% "attractor"
  plotly::layout(p, title = paste(tools::toTitleCase(sysname), "attractor"),
                 scene = list(xaxis = list(title = "x"), yaxis = list(title = "y"),
                              zaxis = list(title = "z"), aspectmode = "data"),
                 legend = list(x = 0, y = 1))
}
