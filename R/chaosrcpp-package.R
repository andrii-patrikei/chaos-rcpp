#' chaosrcpp: chaotic dynamical systems in Rcpp with precision control and RQA
#'
#' Fast compiled integrators and iterators for classic chaotic systems, a
#' floating-point precision switch (fp32, fp64, double-double and an
#' arbitrary-precision Rmpfr reference) for demonstrating sensitive dependence
#' on rounding error,
#' Lyapunov exponents, recurrence plots and Recurrence Quantification Analysis
#' (RQA) measures, Dynamic Mode Decomposition (plain, delay-coordinate and
#' HAVOK), plotting helpers and a small Shiny explorer.
#'
#' Main entry points:
#' \itemize{
#'   \item Flows: [lorenz()], [roessler()], [lyapunov_spectrum()],
#'     [precision_divergence()]
#'   \item Maps: [logistic_bifurcation()], [logistic_lyapunov()],
#'     [logistic_cobweb()], [henon()], [ikeda()], [standard_map()]
#'   \item RQA: [embed_delay()], [recurrence_matrix()], [rqa()],
#'     [rqa_features()]
#'   \item Koopman / DMD: [dmd()], [dmd_features()], [predict.dmd()],
#'     [havok()]
#'   \item Plots: [plot_divergence()], [plot_bifurcation()],
#'     [plot_recurrence()], [plot_attractor()], [plot_attractor_3d()]
#'   \item Benchmark: [benchmark_lorenz()]
#'   \item Shiny: [run_chaos_explorer()]
#' }
#'
#' @docType package
#' @name chaosrcpp-package
#' @aliases chaosrcpp
#' @useDynLib chaosrcpp, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats quantile sd
#' @importFrom graphics plot points lines abline image par layout legend axis box mtext
#' @importFrom stats fitted predict
#' @importFrom grDevices hcl.colors rgb adjustcolor
"_PACKAGE"

# Stop with a helpful message when an optional package is missing.
need_pkg <- function(pkg, purpose) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is needed for %s. Install it with install.packages(\"%s\").",
                 pkg, purpose, pkg), call. = FALSE)
  }
  invisible(TRUE)
}

# Unit roundoff of each compiled precision level.
precision_epsilon <- c(fp32 = 2^-24, fp64 = 2^-53, dd = 2^-104)
