# Dynamic Mode Decomposition (DMD), Hankel (delay-coordinate) DMD and HAVOK.
#
# Everything here is pure R on top of base LAPACK (svd, eigen, qr): the cost
# is one thin SVD of a (features x snapshots) matrix, which is small for the
# low-dimensional systems in this package even with a few hundred delays.

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# States, names, time step and start time from any input the package accepts.
dmd_input <- function(x, dt) {
  t0 <- 0
  if (inherits(x, "chaos_trajectory")) {
    if (is.null(dt)) dt <- attr(x, "dt")
    t0 <- x$t[1]
  } else if (inherits(x, "chaos_map")) {
    if (is.null(dt)) dt <- 1
  }
  X <- as_states(x, m = 1L, tau = 1L)
  if (is.null(colnames(X))) {
    colnames(X) <- if (ncol(X) == 1L) "x" else paste0("x", seq_len(ncol(X)))
  }
  if (is.null(dt)) dt <- 1
  if (anyNA(X)) stop("missing values are not allowed", call. = FALSE)
  list(X = X, dt = dt, t0 = t0)
}

# Stack `delays` time-shifted copies of the state matrix side by side.
# Row i holds x_i, x_{i+1}, ..., x_{i+delays-1} (each a block of ncol(X)
# columns). Row i therefore ends at time index i + delays - 1.
hankel_states <- function(X, delays) {
  delays <- as.integer(delays)
  if (delays < 1L) stop("`delays` must be >= 1", call. = FALSE)
  n <- nrow(X) - delays + 1L
  if (n < 3L) stop(sprintf("series too short for %d delays", delays), call. = FALSE)
  d <- ncol(X)
  H <- matrix(0, n, d * delays)
  for (j in seq_len(delays)) {
    H[, (j - 1L) * d + seq_len(d)] <- X[j:(j + n - 1L), , drop = FALSE]
  }
  H
}

# Least squares for a complex system, Phi %*% b = z.
cls <- function(Phi, z) {
  qr.solve(Phi, as.complex(z), tol = 1e-12)
}

# ---------------------------------------------------------------------------
# dmd()
# ---------------------------------------------------------------------------

#' Dynamic Mode Decomposition
#'
#' Fits the best-fit linear operator \eqn{A} with \eqn{z_{k+1} \approx A z_k}
#' to a sequence of snapshots and returns its eigenvalues and modes (exact DMD
#' of Tu et al. 2014). Each mode \eqn{\phi_j} evolves as \eqn{\lambda_j^k}, so
#' its continuous-time eigenvalue \eqn{\omega_j = \log(\lambda_j)/\Delta t}
#' gives a growth rate \eqn{\mathrm{Re}\,\omega_j} and a frequency
#' \eqn{\mathrm{Im}\,\omega_j / 2\pi}.
#'
#' With `delays > 1` the snapshots are delay vectors (Hankel DMD). This is
#' what makes DMD useful for chaotic systems and for scalar series: a
#' three-variable Lorenz trajectory gives at most three DMD modes, far too
#' few for a linear model, whereas a few dozen delay coordinates approximate
#' the Koopman operator of the attractor (Arbabi and Mezic 2017). A periodic
#' or quasi-periodic signal collapses onto a few isolated eigenvalues on the
#' unit circle. A finely sampled chaotic flow needs many modes, most of them
#' close to the unit circle (the finite-dimensional trace of Koopman's
#' continuous spectrum). A chaotic map, whose correlations decay within a few
#' iterations, needs many modes that sit well inside the circle, and its
#' one-step residual grows with the Lyapunov exponent (see [dmd_features()]).
#'
#' DMD is a diagnostic here, not a forecaster: on the Lorenz attractor a
#' Hankel-DMD model tracks the true trajectory for only about 0.2-0.3 time
#' units, a fraction of the Lyapunov time \eqn{1/\lambda_1 \approx 1.1}.
#'
#' @param x A numeric vector, a matrix or data frame of states (rows are time),
#'   a `chaos_trajectory` or a `chaos_map`.
#' @param rank Truncation rank of the SVD. When `NULL` the smallest rank
#'   capturing a fraction `energy` of the squared singular values is used.
#' @param delays Number of delay copies stacked into each snapshot; `1` is
#'   plain DMD.
#' @param dt Sampling interval. Taken from a `chaos_trajectory`, `1` for
#'   maps and plain vectors.
#' @param energy Energy threshold used when `rank = NULL`.
#' @param center Subtract the column means before fitting. Mean-subtracted
#'   DMD at full rank reduces to a temporal DFT (Chen et al. 2012), so leave
#'   this `FALSE` unless you truncate.
#' @return An object of class `dmd`: a list with `eigenvalues`
#'   (discrete-time \eqn{\lambda}), `omega`, `growth`, `frequency`, `modes`
#'   (columns, sorted by energy), `amplitudes` (at the first snapshot),
#'   `amplitudes_end` (at the last snapshot, used by [predict.dmd()]),
#'   `energy` (normalised mode energies), `Atilde` (the reduced operator),
#'   `singular_values`, `one_step_nrmse`, and the settings.
#' @references
#' Schmid, P. J. (2010). Dynamic mode decomposition of numerical and
#' experimental data. Journal of Fluid Mechanics, 656, 5-28.
#'
#' Tu, J. H., Rowley, C. W., Luchtenburg, D. M., Brunton, S. L., & Kutz, J. N.
#' (2014). On dynamic mode decomposition: theory and applications. Journal of
#' Computational Dynamics, 1(2), 391-421.
#'
#' Arbabi, H., & Mezic, I. (2017). Ergodic theory, dynamic mode decomposition,
#' and computation of spectral properties of the Koopman operator. SIAM
#' Journal on Applied Dynamical Systems, 16(4), 2096-2126.
#' @examples
#' # two damped oscillations are recovered exactly
#' t <- seq(0, 20, by = 0.05)
#' x <- exp(-0.1 * t) * sin(2 * pi * 0.5 * t) + 0.5 * sin(2 * pi * 1.3 * t)
#' d <- dmd(x, rank = 4, delays = 20, dt = 0.05)
#' d
#' plot(d)
#'
#' # a chaotic attractor: eigenvalues crowd onto the unit circle
#' tr <- lorenz(n = 8000, dt = 0.01)
#' plot(dmd(tr, delays = 60, rank = 40))
#' @seealso [dmd_features()], [predict.dmd()], [havok()]
#' @export
dmd <- function(x, rank = NULL, delays = 1L, dt = NULL, energy = 0.999,
                center = FALSE) {
  inp <- dmd_input(x, dt)
  X <- inp$X
  mu <- if (center) colMeans(X) else rep(0, ncol(X))
  Xc <- sweep(X, 2L, mu)
  H <- hankel_states(Xc, delays)
  Z <- t(H)                                   # features x snapshots
  X1 <- Z[, -ncol(Z), drop = FALSE]
  X2 <- Z[, -1L, drop = FALSE]

  s <- svd(X1)
  sv <- s$d
  numrank <- sum(sv > max(sv) * max(dim(X1)) * .Machine$double.eps)
  if (numrank == 0L) stop("the snapshot matrix is zero", call. = FALSE)
  if (is.null(rank)) {
    ce <- cumsum(sv^2) / sum(sv^2)
    rank <- which(ce >= energy)[1L]
  }
  rank <- as.integer(min(rank, numrank))
  U <- s$u[, seq_len(rank), drop = FALSE]
  V <- s$v[, seq_len(rank), drop = FALSE]
  sinv <- 1 / sv[seq_len(rank)]

  X2V <- X2 %*% V
  Atilde <- crossprod(U, X2V) * rep(sinv, each = rank)   # U' X2 V S^-1
  e <- eigen(Atilde)
  lambda <- e$values
  W <- e$vectors
  Phi <- X2V %*% (sinv * W)                                # exact DMD modes

  b0 <- cls(Phi, Z[, 1L])
  bN <- cls(Phi, Z[, ncol(Z)])
  en <- Mod(b0)^2 * colSums(Mod(Phi)^2)
  ord <- order(en, decreasing = TRUE)

  omega <- log(as.complex(lambda)) / inp$dt
  resid <- X2 - U %*% (Atilde %*% crossprod(U, X1))
  nrmse <- sqrt(sum(resid^2) / sum(X2^2))

  out <- list(
    eigenvalues = lambda[ord],
    omega = omega[ord],
    growth = Re(omega[ord]),
    frequency = Im(omega[ord]) / (2 * pi),
    modes = Phi[, ord, drop = FALSE],
    amplitudes = b0[ord],
    amplitudes_end = bN[ord],
    energy = en[ord] / sum(en),
    Atilde = Atilde,
    singular_values = sv,
    one_step_nrmse = nrmse,
    rank = rank, delays = as.integer(delays), dt = inp$dt, t0 = inp$t0,
    n_vars = ncol(X), var_names = colnames(X), center = mu,
    n_snapshots = ncol(Z), n_time = nrow(X)
  )
  class(out) <- "dmd"
  out
}

#' @export
print.dmd <- function(x, n = 6L, digits = 4, ...) {
  cat(sprintf("<dmd> rank %d, %d delay%s, %d variable%s, %d snapshots, dt = %s\n",
              x$rank, x$delays, if (x$delays == 1L) "" else "s", x$n_vars,
              if (x$n_vars == 1L) "" else "s", x$n_snapshots, format(x$dt)))
  cat(sprintf("  one-step relative error %.3g; spectral radius %.4f\n",
              x$one_step_nrmse, max(Mod(x$eigenvalues))))
  k <- seq_len(min(n, x$rank))
  tab <- data.frame(`|lambda|` = Mod(x$eigenvalues[k]), growth = x$growth[k],
                    frequency = x$frequency[k], energy = x$energy[k],
                    check.names = FALSE)
  print(signif(tab, digits))
  if (x$rank > n) cat(sprintf("  ... %d more modes\n", x$rank - n))
  invisible(x)
}

# Map snapshot indices to the original variables. Snapshot j ends at time
# index j + delays - 1, so the last block of each delay vector is "now".
dmd_snapshots_to_series <- function(object, Zhat, idx) {
  d <- object$n_vars
  last <- (object$delays - 1L) * d + seq_len(d)
  Y <- t(Re(Zhat[last, , drop = FALSE]))
  Y <- sweep(Y, 2L, object$center, "+")
  colnames(Y) <- object$var_names
  time_index <- idx + object$delays - 1L
  data.frame(t = object$t0 + (time_index - 1) * object$dt, Y)
}

#' Reconstruct or forecast with a DMD model
#'
#' `fitted()` rebuilds the training series from the first snapshot,
#' \eqn{z_k = \sum_j b_j \lambda_j^{k} \phi_j}. `predict()` continues it
#' `n_ahead` steps past the last observation, starting from amplitudes fitted
#' to the final snapshot. For chaotic input expect useful forecasts only for
#' a fraction of a Lyapunov time: the model is linear and, with all
#' \eqn{|\lambda| \le 1}, it can only reproduce the quasi-periodic part of the
#' motion.
#'
#' @param object A `dmd` object.
#' @param n_ahead Number of steps to forecast.
#' @param ... Unused.
#' @return A data frame with `t` and one column per original variable.
#' @examples
#' t <- seq(0, 20, by = 0.05)
#' x <- sin(2 * pi * 0.4 * t) + 0.3 * sin(2 * pi * 1.1 * t)
#' d <- dmd(x, rank = 4, delays = 10, dt = 0.05)
#' head(predict(d, n_ahead = 100))
#' @export
predict.dmd <- function(object, n_ahead = 100L, ...) {
  k <- seq_len(as.integer(n_ahead))
  Zhat <- object$modes %*% (object$amplitudes_end * outer(object$eigenvalues, k, `^`))
  dmd_snapshots_to_series(object, Zhat, object$n_snapshots + k)
}

#' @rdname predict.dmd
#' @export
fitted.dmd <- function(object, ...) {
  k <- seq_len(object$n_snapshots) - 1L
  Zhat <- object$modes %*% (object$amplitudes * outer(object$eigenvalues, k, `^`))
  dmd_snapshots_to_series(object, Zhat, k + 1L)
}

#' Plot a DMD spectrum
#'
#' Left: discrete-time eigenvalues with the unit circle, point area
#' proportional to mode energy. Right: mode energy against frequency.
#'
#' @param x A `dmd` object.
#' @param col Point colour.
#' @param ... Passed to the first [graphics::plot()] call.
#' @return `x`, invisibly.
#' @export
plot.dmd <- function(x, col = "#1f2a44", ...) {
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
  on.exit(par(op))
  lam <- x$eigenvalues
  cex <- 0.6 + 3 * sqrt(x$energy / max(x$energy))
  lim <- max(1.05, max(Mod(lam)) * 1.05)
  th <- seq(0, 2 * pi, length.out = 361)
  plot(cos(th), sin(th), type = "l", col = "grey60", asp = 1,
       xlim = c(-lim, lim), ylim = c(-lim, lim),
       xlab = "Re(lambda)", ylab = "Im(lambda)", main = "DMD eigenvalues", ...)
  abline(h = 0, v = 0, col = "grey85")
  points(Re(lam), Im(lam), pch = 16, cex = cex, col = adjustcolor(col, 0.7))
  keep <- x$frequency >= 0
  plot(x$frequency[keep], x$energy[keep], type = "h", lwd = 2, col = col,
       xlab = if (x$dt == 1) "frequency (cycles per step)" else "frequency",
       ylab = "energy", main = "Mode spectrum")
  points(x$frequency[keep], x$energy[keep], pch = 16, col = col)
  invisible(x)
}

# ---------------------------------------------------------------------------
# dmd_features()
# ---------------------------------------------------------------------------

dmd_feature_names <- c("dmd_rank", "dmd_dom_freq", "dmd_dom_growth", "dmd_max_growth",
                       "dmd_spectral_radius", "dmd_unit_circle_frac", "dmd_top_energy",
                       "dmd_spectral_entropy", "dmd_one_step_nrmse")

#' DMD summary measures as a one-row data frame
#'
#' The DMD counterpart of [rqa_features()]: a fixed set of scalars per
#' signal or window, ready for a feature table.
#'
#' * `dmd_rank`: rank used (with `rank = NULL`, how many modes are needed to
#'   reach the energy threshold - a linear-algebra complexity measure).
#' * `dmd_dom_freq`, `dmd_dom_growth`: frequency and growth rate of the most
#'   energetic oscillating mode.
#' * `dmd_max_growth`: largest growth rate (positive means unstable modes).
#' * `dmd_spectral_radius`: largest \eqn{|\lambda|}.
#' * `dmd_unit_circle_frac`: fraction of eigenvalues within `tol` of the unit
#'   circle: 1 for periodic and quasi-periodic signals, high for finely
#'   sampled flows, near 0 for chaotic maps.
#' * `dmd_top_energy`: energy share of the leading mode.
#' * `dmd_spectral_entropy`: Shannon entropy of the mode energies divided by
#'   `log(rank)` (0 = one mode, 1 = energy spread evenly).
#' * `dmd_one_step_nrmse`: relative one-step residual of the linear model; low
#'   for regular dynamics, higher for chaotic or noisy ones. Over a logistic
#'   map sweep \eqn{r \in [3.4, 4]} (16 delays, centred) its Spearman
#'   correlation with the Lyapunov exponent is about 0.96, and it drops back
#'   inside the period-3 window.
#'
#' Frequencies and growth rates are per unit of `dt`, so compare features
#' only between signals sampled at the same rate.
#'
#' @inheritParams dmd
#' @param tol Tolerance for `dmd_unit_circle_frac`.
#' @param ... Passed to [dmd()].
#' @return A one-row data frame.
#' @examples
#' dmd_features(sin(seq(0, 30 * pi, length.out = 1000)), delays = 10, rank = 4)
#' @export
dmd_features <- function(x, ..., tol = 0.01) {
  d <- dmd(x, ...)
  osc <- which(abs(d$frequency) > 1e-8)
  dom <- if (length(osc)) osc[1L] else 1L
  p <- d$energy[d$energy > 0]
  data.frame(
    dmd_rank = d$rank,
    dmd_dom_freq = abs(d$frequency[dom]),
    dmd_dom_growth = d$growth[dom],
    dmd_max_growth = max(d$growth),
    dmd_spectral_radius = max(Mod(d$eigenvalues)),
    dmd_unit_circle_frac = mean(abs(Mod(d$eigenvalues) - 1) < tol),
    dmd_top_energy = d$energy[1L],
    dmd_spectral_entropy = if (d$rank > 1L) -sum(p * log(p)) / log(d$rank) else 0,
    dmd_one_step_nrmse = d$one_step_nrmse
  )
}

# ---------------------------------------------------------------------------
# havok()
# ---------------------------------------------------------------------------

#' Hankel Alternative View Of Koopman (HAVOK)
#'
#' Brunton et al. (2017) showed that a chaotic scalar series, viewed through
#' enough delay coordinates, is well described by a *linear* system driven by
#' a single intermittent forcing term:
#' \deqn{\dot v = A v + B v_r,}
#' where \eqn{v_1, \dots, v_r} are the leading right singular vectors of the
#' Hankel matrix (eigen-time-delay coordinates) and the last one, \eqn{v_r},
#' plays the role of the forcing. For the Lorenz system the fitted linear
#' part is almost exactly skew-symmetric (a chain of oscillators) and the
#' forcing is mostly quiet, bursting in the saddle region where lobe
#' switching happens. The association is one-sided. Starting on the attractor
#' with the paper's settings (`dt = 0.001`, 100 delays, rank 15), about 93% of
#' forcing bursts above 2 sd fall within 0.3 time units of a lobe switch, but
#' only about a third of switches come with such a burst. A burst is a
#' reliable warning; its absence is not an all-clear. Discard the initial
#' transient first, or it produces spurious bursts.
#'
#' Derivatives of the coordinates are estimated by fourth-order central
#' differences and the model is fitted by ordinary least squares (the
#' original paper uses sparse regression; the fitted `A` is close to
#' skew-symmetric and tridiagonal either way).
#'
#' @param x A numeric vector or a `chaos_trajectory` (then `var` selects the
#'   measured coordinate).
#' @param delays Number of delays (rows of the Hankel matrix). Brunton et al.
#'   use a window of about 0.1 time units for Lorenz.
#' @param rank Number of delay coordinates kept; the last is the forcing.
#' @param dt Sampling interval (taken from a `chaos_trajectory`).
#' @param var Column used when `x` is a `chaos_trajectory`.
#' @return An object of class `havok`: `A` (\eqn{(r-1) \times (r-1)}), `B`,
#'   `forcing` (\eqn{v_r}), `coords` (the \eqn{v} matrix), `t`, `x`
#'   (aligned measured series), `singular_values`, `fit_r2` and settings.
#' @references Brunton, S. L., Brunton, B. W., Proctor, J. L., Kaiser, E., &
#'   Kutz, J. N. (2017). Chaos as an intermittently forced linear system.
#'   Nature Communications, 8, 19.
#' @examples
#' tr <- lorenz(n = 40000, dt = 0.002)
#' h <- havok(tr, delays = 50, rank = 11)
#' plot(h)
#' @seealso [dmd()]
#' @export
havok <- function(x, delays = 100L, rank = 15L, dt = NULL, var = "x") {
  t0 <- 0
  if (inherits(x, "chaos_trajectory")) {
    if (is.null(dt)) dt <- attr(x, "dt")
    t0 <- x$t[1]
    x <- x[[var]]
  }
  if (is.null(dt)) dt <- 1
  x <- as.numeric(x)
  rank <- as.integer(rank)
  if (rank < 2L || rank > delays) stop("need 2 <= rank <= delays", call. = FALSE)

  H <- hankel_states(matrix(x, ncol = 1L), delays)        # snapshots x delays
  s <- svd(H, nu = rank, nv = 0L)
  Vr <- s$u                                               # time series of coordinates
  n <- nrow(Vr)
  i <- 3:(n - 2L)
  dV <- (-Vr[i + 2L, ] + 8 * Vr[i + 1L, ] - 8 * Vr[i - 1L, ] + Vr[i - 2L, ]) / (12 * dt)
  Vi <- Vr[i, , drop = FALSE]
  Xi <- qr.solve(Vi, dV[, seq_len(rank - 1L), drop = FALSE])
  A <- t(Xi[seq_len(rank - 1L), , drop = FALSE])
  B <- Xi[rank, ]
  pred <- Vi %*% Xi
  target <- dV[, seq_len(rank - 1L), drop = FALSE]
  r2 <- 1 - sum((target - pred)^2) / sum(sweep(target, 2L, colMeans(target))^2)

  # align: snapshot j of the Hankel matrix ends at time index j + delays - 1
  idx <- seq_len(n) + delays - 1L
  out <- list(A = A, B = B, forcing = Vr[, rank], coords = Vr,
              t = t0 + (idx - 1) * dt, x = x[idx],
              singular_values = s$d, fit_r2 = r2,
              delays = as.integer(delays), rank = rank, dt = dt)
  class(out) <- "havok"
  out
}

#' @export
print.havok <- function(x, ...) {
  S <- x$A + t(x$A)
  cat(sprintf("<havok> %d delays, rank %d (forcing = v%d), dt = %s\n",
              x$delays, x$rank, x$rank, format(x$dt)))
  cat(sprintf("  derivative fit R^2 = %.4f; |A + A'| / |A| = %.3f (0 = skew-symmetric)\n",
              x$fit_r2, sqrt(sum(S^2)) / sqrt(sum(x$A^2))))
  cat(sprintf("  forcing active (|v_r| > 2 sd) %.1f%% of the time\n",
              100 * mean(abs(x$forcing) > 2 * stats::sd(x$forcing))))
  invisible(x)
}

#' @rdname havok
#' @param threshold Forcing magnitude, in standard deviations, above which the
#'   forcing is drawn as active.
#' @param col Colour of active-forcing points.
#' @param ... Unused.
#' @export
plot.havok <- function(x, threshold = 2, col = "#c0392b", ...) {
  op <- par(mfrow = c(2, 1), mar = c(2.5, 4, 2, 1), oma = c(1.5, 0, 0, 0))
  on.exit(par(op))
  active <- abs(x$forcing) > threshold * stats::sd(x$forcing)
  plot(x$t, x$x, type = "l", col = "grey40", xlab = "", ylab = "measured",
       main = "Measured series, active forcing in colour")
  points(x$t[active], x$x[active], pch = 16, cex = 0.3, col = col)
  plot(x$t, x$forcing, type = "l", col = "grey40", xlab = "",
       ylab = bquote(v[.(x$rank)]), main = "HAVOK forcing")
  abline(h = c(-1, 1) * threshold * stats::sd(x$forcing), lty = 2, col = col)
  mtext("t", side = 1, outer = TRUE, line = 0.3)
  invisible(x)
}
