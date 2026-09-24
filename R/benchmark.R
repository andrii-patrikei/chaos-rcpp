#' Benchmark the compiled integrator against deSolve and pure R
#'
#' Runs the same RK4 integration of the Lorenz system (same `dt`, same
#' number of steps) with four implementations and times each:
#' \itemize{
#'   \item `chaosrcpp fp64` and `chaosrcpp fp32`: the templated C++ kernel
#'   \item `deSolve rk4`: [deSolve::ode()] with `method = "rk4"` and an R
#'     right-hand-side function (skipped when deSolve is not installed)
#'   \item `pure R`: an explicit RK4 loop written in R
#' }
#' All four implement the identical scheme, so the trajectories agree to
#' rounding over a short horizon; the column `max_abs_diff` reports the
#' largest deviation from `chaosrcpp fp64` over the first `compare_steps`
#' steps as a correctness check.
#'
#' @param n Number of RK4 steps.
#' @param dt Time step.
#' @param x0 Initial state.
#' @param compare_steps Steps over which trajectories are compared.
#' @param pure_r Include the pure R loop (slow for large `n`).
#' @return Data frame with columns `method`, `seconds`, `steps_per_second`,
#'   `speedup_vs_deSolve` and `max_abs_diff`, class `chaos_benchmark`.
#' @examples
#' \donttest{
#' benchmark_lorenz(n = 20000)
#' }
#' @export
benchmark_lorenz <- function(n = 100000L, dt = 0.01, x0 = c(1, 1, 1),
                             compare_steps = 500L, pure_r = TRUE) {
  n <- as.integer(n)
  sigma <- 10; rho <- 28; beta <- 8 / 3
  timings <- list()
  results <- list()

  time_it <- function(expr) {
    st <- system.time(val <- expr)
    list(seconds = unname(st[["elapsed"]]), value = val)
  }

  r <- time_it(lorenz(x0, sigma, rho, beta, dt, n, precision = "fp64"))
  timings[["chaosrcpp fp64"]] <- r$seconds
  results[["chaosrcpp fp64"]] <- as.matrix(as.data.frame(unclass(r$value))[, c("x", "y", "z")])

  r <- time_it(lorenz(x0, sigma, rho, beta, dt, n, precision = "fp32"))
  timings[["chaosrcpp fp32"]] <- r$seconds
  results[["chaosrcpp fp32"]] <- as.matrix(as.data.frame(unclass(r$value))[, c("x", "y", "z")])

  if (requireNamespace("deSolve", quietly = TRUE)) {
    rhs <- function(t, state, p) {
      list(c(p$sigma * (state[2] - state[1]),
             state[1] * (p$rho - state[3]) - state[2],
             state[1] * state[2] - p$beta * state[3]))
    }
    times <- seq(0, by = dt, length.out = n + 1)
    r <- time_it(deSolve::ode(y = x0, times = times, func = rhs,
                              parms = list(sigma = sigma, rho = rho, beta = beta),
                              method = "rk4"))
    timings[["deSolve rk4"]] <- r$seconds
    results[["deSolve rk4"]] <- unname(as.matrix(r$value[, 2:4]))
  }

  if (pure_r) {
    r <- time_it(lorenz_rk4_pure_r(x0, sigma, rho, beta, dt, n))
    timings[["pure R"]] <- r$seconds
    results[["pure R"]] <- r$value
  }

  ref <- results[["chaosrcpp fp64"]][seq_len(compare_steps + 1L), ]
  diffs <- vapply(names(results), function(k) {
    max(abs(results[[k]][seq_len(compare_steps + 1L), ] - ref))
  }, numeric(1))

  secs <- unlist(timings)
  out <- data.frame(method = names(secs), seconds = unname(secs),
                    steps_per_second = n / unname(secs),
                    speedup_vs_deSolve = if ("deSolve rk4" %in% names(secs)) unname(secs[["deSolve rk4"]] / secs) else NA_real_,
                    max_abs_diff = unname(diffs[names(secs)]),
                    stringsAsFactors = FALSE)
  attr(out, "n") <- n
  attr(out, "dt") <- dt
  class(out) <- c("chaos_benchmark", "data.frame")
  out
}

# Reference implementation: RK4 for Lorenz in plain R.
lorenz_rk4_pure_r <- function(x0, sigma, rho, beta, dt, n) {
  f <- function(v) c(sigma * (v[2] - v[1]), v[1] * (rho - v[3]) - v[2], v[1] * v[2] - beta * v[3])
  out <- matrix(0, n + 1, 3)
  x <- as.numeric(x0)
  out[1, ] <- x
  for (k in seq_len(n)) {
    k1 <- f(x)
    k2 <- f(x + 0.5 * dt * k1)
    k3 <- f(x + 0.5 * dt * k2)
    k4 <- f(x + dt * k3)
    x <- x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    out[k + 1, ] <- x
  }
  out
}

#' @export
print.chaos_benchmark <- function(x, ...) {
  cat(sprintf("Lorenz RK4, %d steps, dt = %g\n", attr(x, "n"), attr(x, "dt")))
  df <- as.data.frame(unclass(x))
  df$seconds <- signif(df$seconds, 3)
  df$steps_per_second <- format(round(df$steps_per_second), big.mark = ",")
  df$speedup_vs_deSolve <- ifelse(is.na(df$speedup_vs_deSolve), "", sprintf("%.0fx", df$speedup_vs_deSolve))
  df$max_abs_diff <- format(df$max_abs_diff, digits = 2)
  print(df, row.names = FALSE, ...)
  invisible(x)
}
