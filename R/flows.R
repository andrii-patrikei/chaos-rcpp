#' Integrate the Lorenz system
#'
#' Classical fourth-order Runge-Kutta integration of the Lorenz equations
#' \deqn{\dot x = \sigma (y - x), \quad \dot y = x (\rho - z) - y, \quad
#'       \dot z = x y - \beta z}
#' in a chosen floating-point precision.
#'
#' `"fp32"`, `"fp64"` and `"dd"` run the same templated C++ code instantiated
#' with `float`, `double` and a double-double type (two doubles carrying about
#' 106 bits, roughly 32 significant digits, built from error-free
#' transformations). `"mpfr"` runs the same RK4 formulas in arbitrary
#' precision through the Rmpfr package; it is the gold-standard reference but
#' runs as an R loop at roughly a hundred steps per second, so use it for
#' validation and short runs and `"dd"` for long high-precision runs. Because
#' all levels share the algorithm and the step size, differences between them
#' isolate the effect of rounding error. Note that even the high-precision runs
#' carry the truncation error of RK4 itself, so they are a reference for the
#' *rounding* question, not the exact solution of the ODE.
#'
#' @param x0 Initial state, numeric vector of length 3.
#' @param sigma,rho,beta Lorenz parameters. The classic chaotic set is
#'   `sigma = 10, rho = 28, beta = 8/3`.
#' @param dt Time step.
#' @param n Number of RK4 steps.
#' @param precision One of `"fp64"`, `"fp32"`, `"dd"` or `"mpfr"`.
#' @param thin Store every `thin`-th step (useful for very long runs).
#' @param mpfr_bits Mantissa bits for the `"mpfr"` precision.
#' @return A data frame with columns `t`, `x`, `y`, `z` and class
#'   `chaos_trajectory`. Attributes `system`, `precision`, `dt` and `params`
#'   describe the run.
#' @examples
#' tr <- lorenz(n = 5000)
#' plot(tr)
#' @seealso [roessler()], [precision_divergence()], [lyapunov_spectrum()]
#' @export
lorenz <- function(x0 = c(1, 1, 1), sigma = 10, rho = 28, beta = 8 / 3,
                   dt = 0.01, n = 10000L,
                   precision = c("fp64", "fp32", "dd", "mpfr"),
                   thin = 1L, mpfr_bits = 256L) {
  precision <- match.arg(precision)
  n <- as.integer(n)
  thin <- as.integer(thin)
  params <- list(sigma = sigma, rho = rho, beta = beta)

  if (precision == "mpfr") {
    rhs <- function(x, p) {
      c(p$sigma * (x[2] - x[1]),
        x[1] * (p$rho - x[3]) - x[2],
        x[1] * x[2] - p$beta * x[3])
    }
    m <- flow_mpfr(rhs, x0, params, dt, n, thin, mpfr_bits)
  } else {
    m <- lorenz_cpp(as.numeric(x0), sigma, rho, beta, dt, n, thin,
                    precision_code(precision))
  }
  new_trajectory(m, dt * thin, "lorenz", precision, params)
}

#' Integrate the Roessler system
#'
#' RK4 integration of
#' \deqn{\dot x = -y - z, \quad \dot y = x + a y, \quad \dot z = b + z (x - c)}
#' with the same precision options as [lorenz()].
#'
#' @inheritParams lorenz
#' @param a,b,c Roessler parameters; `a = 0.2, b = 0.2, c = 5.7` is the
#'   classic chaotic set.
#' @return A `chaos_trajectory` data frame, see [lorenz()].
#' @examples
#' tr <- roessler(n = 20000)
#' plot(tr, vars = c("x", "y"))
#' @export
roessler <- function(x0 = c(1, 1, 0), a = 0.2, b = 0.2, c = 5.7,
                     dt = 0.02, n = 20000L,
                     precision = c("fp64", "fp32", "dd", "mpfr"),
                     thin = 1L, mpfr_bits = 256L) {
  precision <- match.arg(precision)
  n <- as.integer(n)
  thin <- as.integer(thin)
  params <- list(a = a, b = b, c = c)

  if (precision == "mpfr") {
    rhs <- function(x, p) {
      c(-x[2] - x[3],
        x[1] + p$a * x[2],
        p$b + x[3] * (x[1] - p$c))
    }
    m <- flow_mpfr(rhs, x0, params, dt, n, thin, mpfr_bits)
  } else {
    m <- roessler_cpp(as.numeric(x0), a, b, c, dt, n, thin,
                      precision_code(precision))
  }
  new_trajectory(m, dt * thin, "roessler", precision, params)
}

precision_code <- function(precision) {
  c(fp64 = 0L, fp32 = 1L, dd = 2L)[[precision]]
}

# RK4 in arbitrary precision with Rmpfr. Same formulas as the C++ kernel.
# Returns the double-rounded states as a matrix, with the exact mpfrMatrix
# attached as attribute "mpfr".
flow_mpfr <- function(rhs, x0, params, dt, n, thin, bits) {
  need_pkg("Rmpfr", "the \"mpfr\" precision")
  bits <- as.integer(bits)
  x <- Rmpfr::mpfr(x0, bits)
  p <- lapply(params, Rmpfr::mpfr, precBits = bits)
  h <- Rmpfr::mpfr(dt, bits)
  half_h <- h / 2
  sixth_h <- h / 6

  n_out <- n %/% thin + 1L
  kept <- vector("list", n_out)
  kept[[1L]] <- x
  row <- 2L
  for (k in seq_len(n)) {
    k1 <- rhs(x, p)
    k2 <- rhs(x + half_h * k1, p)
    k3 <- rhs(x + half_h * k2, p)
    k4 <- rhs(x + h * k3, p)
    x <- x + sixth_h * (k1 + 2 * k2 + 2 * k3 + k4)
    if (k %% thin == 0L && row <= n_out) {
      kept[[row]] <- x
      row <- row + 1L
    }
  }
  exact <- t(Rmpfr::mpfr2array(do.call(c, kept), dim = c(length(x0), n_out)))
  out <- matrix(as.numeric(exact), n_out, length(x0))
  attr(out, "mpfr") <- exact
  out
}

new_trajectory <- function(m, dt_out, system, precision, params) {
  D <- 3L
  df <- data.frame(t = seq(0, by = dt_out, length.out = nrow(m)),
                   x = m[, 1], y = m[, 2], z = m[, 3])
  attr(df, "system") <- system
  attr(df, "precision") <- precision
  attr(df, "dt") <- dt_out
  attr(df, "params") <- params
  if (precision == "dd") {
    # low-order parts of the double-double states; hi + lo is the exact value
    attr(df, "lo") <- unname(m[, D + seq_len(D), drop = FALSE])
  }
  if (precision == "mpfr") attr(df, "mpfr") <- attr(m, "mpfr")
  class(df) <- c("chaos_trajectory", "data.frame")
  df
}

#' Exact state values of a high-precision trajectory
#'
#' `"fp32"` and `"fp64"` trajectories are exactly what the data frame holds.
#' `"dd"` and `"mpfr"` trajectories carry more digits than a double, so the
#' data frame columns are rounded and the full values live in attributes.
#' This helper returns the states as an Rmpfr matrix (requires Rmpfr) so
#' that differences between precisions can be computed without losing them
#' to the rounding of the output itself.
#'
#' @param x A `chaos_trajectory`.
#' @param bits Precision of the returned mpfr numbers (ignored for `"mpfr"`
#'   trajectories, which keep their own).
#' @return An `mpfrMatrix` with one row per state and columns x, y, z.
#' @examples
#' if (requireNamespace("Rmpfr", quietly = TRUE)) {
#'   tr <- lorenz(n = 5, precision = "dd")
#'   exact_states(tr)
#' }
#' @export
exact_states <- function(x, bits = 128L) {
  need_pkg("Rmpfr", "exact_states()")
  if (!is.null(attr(x, "mpfr"))) return(attr(x, "mpfr"))
  M <- Rmpfr::mpfr(as.matrix(as.data.frame(unclass(x))[, c("x", "y", "z")]), bits)
  lo <- attr(x, "lo")
  if (!is.null(lo)) M <- M + Rmpfr::mpfr(lo, bits)
  M
}

#' @export
print.chaos_trajectory <- function(x, ...) {
  p <- attr(x, "params")
  cat(sprintf("<chaos_trajectory> %s system, precision %s, %d states, dt = %g\n",
              attr(x, "system"), attr(x, "precision"), nrow(x), attr(x, "dt")))
  cat("  parameters:", paste(names(p), signif(unlist(p), 6), sep = " = ", collapse = ", "), "\n")
  print(utils::head(as.data.frame(unclass(x)), 6), ...)
  if (nrow(x) > 6) cat(sprintf("  ... %d more rows\n", nrow(x) - 6))
  invisible(x)
}

#' Lyapunov spectrum of a chaotic flow
#'
#' Estimates all Lyapunov exponents with the standard method of Benettin and
#' co-workers: the state and a full set of tangent vectors are integrated
#' together with RK4 and the tangent vectors are re-orthonormalised by
#' Gram-Schmidt after every step. The average logarithmic growth rate of the
#' `j`-th orthogonalised vector is the `j`-th exponent.
#'
#' For the classic Lorenz parameters the spectrum is approximately
#' `(0.906, 0, -14.57)`, and the exponents must sum to the trace of the
#' Jacobian, `-(sigma + 1 + beta)`, which is a handy self-check.
#'
#' @param system `"lorenz"` or `"roessler"`.
#' @param params Named list of parameters; defaults to the classic chaotic set.
#' @param x0 Initial state (length 3). Defaults to a point near the attractor.
#' @param dt Time step.
#' @param n Number of steps used for averaging.
#' @param n_transient Steps discarded before the tangent vectors start.
#' @return Numeric vector of exponents, sorted from largest to smallest, with
#'   attributes `sum` and `kaplan_yorke` (the Kaplan-Yorke fractal dimension).
#' @examples
#' ly <- lyapunov_spectrum("lorenz", n = 50000)
#' ly
#' attr(ly, "kaplan_yorke")
#' @export
lyapunov_spectrum <- function(system = c("lorenz", "roessler"), params = NULL,
                              x0 = NULL, dt = 0.01, n = 100000L,
                              n_transient = 5000L) {
  system <- match.arg(system)
  n <- as.integer(n)
  n_transient <- as.integer(n_transient)
  if (system == "lorenz") {
    p <- utils::modifyList(list(sigma = 10, rho = 28, beta = 8 / 3), params %||% list())
    if (is.null(x0)) x0 <- c(1, 1, 1)
    ly <- lorenz_lyapunov_cpp(as.numeric(x0), p$sigma, p$rho, p$beta, dt, n, n_transient)
  } else {
    p <- utils::modifyList(list(a = 0.2, b = 0.2, c = 5.7), params %||% list())
    if (is.null(x0)) x0 <- c(1, 1, 0)
    ly <- roessler_lyapunov_cpp(as.numeric(x0), p$a, p$b, p$c, dt, n, n_transient)
  }
  ly <- sort(ly, decreasing = TRUE)
  names(ly) <- paste0("lambda", seq_along(ly))
  attr(ly, "sum") <- sum(ly)
  attr(ly, "kaplan_yorke") <- kaplan_yorke_dimension(ly)
  attr(ly, "system") <- system
  attr(ly, "params") <- p
  class(ly) <- c("lyapunov_spectrum", "numeric")
  ly
}

#' Kaplan-Yorke dimension from a Lyapunov spectrum
#'
#' \eqn{D_{KY} = j + (\lambda_1 + \dots + \lambda_j) / |\lambda_{j+1}|}, with
#' `j` the largest index for which the partial sum is non-negative.
#'
#' @param lambda Lyapunov exponents sorted in decreasing order.
#' @return A single number (NA when the spectrum has no negative exponent).
#' @export
kaplan_yorke_dimension <- function(lambda) {
  lambda <- sort(as.numeric(lambda), decreasing = TRUE)
  cs <- cumsum(lambda)
  j <- max(which(cs >= 0), 0L)
  if (j == 0L) return(0)
  if (j >= length(lambda)) return(NA_real_)
  j + cs[j] / abs(lambda[j + 1])
}

#' @export
print.lyapunov_spectrum <- function(x, digits = 4, ...) {
  cat(sprintf("Lyapunov spectrum (%s):\n", attr(x, "system")))
  print(round(as.numeric(x), digits))
  cat(sprintf("  sum = %.*f   Kaplan-Yorke dimension = %.*f\n",
              digits, attr(x, "sum"), digits, attr(x, "kaplan_yorke")))
  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
