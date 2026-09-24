#' Launch the interactive chaos explorer
#'
#' A small Shiny app with three tabs: the Lorenz attractor in 3D with live
#' parameter sliders and an fp32-versus-fp64 divergence view, the logistic
#' map with a movable `r` slider linked to the bifurcation diagram, cobweb and
#' orbit, and the Chirikov standard map with adjustable kick strength.
#'
#' @param ... Passed on to [shiny::runApp()].
#' @return Called for its side effect.
#' @examples
#' if (interactive()) run_chaos_explorer()
#' @export
run_chaos_explorer <- function(...) {
  need_pkg("shiny", "run_chaos_explorer()")
  need_pkg("plotly", "the 3D view in run_chaos_explorer()")
  dir <- system.file("shiny", "chaos-explorer", package = "chaosrcpp")
  if (dir == "") stop("Shiny app directory not found; reinstall chaosrcpp.", call. = FALSE)
  shiny::runApp(dir, ...)
}
