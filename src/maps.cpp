// Discrete-time chaotic maps. These loops are trivially cheap per iterate,
// which is exactly where an interpreted language hurts most: a bifurcation
// diagram needs millions of iterations and pure R takes seconds to minutes
// for what compiled code does in milliseconds.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

namespace {
const double TWO_PI = 6.283185307179586476925286766559;

inline double wrap_2pi(double v) {
  v = std::fmod(v, TWO_PI);
  if (v < 0.0) v += TWO_PI;
  return v;
}
}  // namespace

// ---------------------------------------------------------------------------
// Logistic map  x_{n+1} = r x_n (1 - x_n)
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector logistic_orbit_cpp(double x0, double r, int n) {
  if (n < 0) stop("n must be non-negative");
  NumericVector out(n + 1);
  double x = x0;
  out[0] = x;
  for (int k = 1; k <= n; ++k) {
    x = r * x * (1.0 - x);
    out[k] = x;
  }
  return out;
}

// For each r: discard n_transient iterates, then keep n_keep of them.
// Returns a list with two long vectors (r repeated, x) ready for plotting.
// [[Rcpp::export]]
List logistic_bifurcation_cpp(NumericVector r, double x0, int n_transient,
                              int n_keep) {
  const R_xlen_t nr = r.size();
  const R_xlen_t total = nr * static_cast<R_xlen_t>(n_keep);
  NumericVector rr(total), xx(total);
  R_xlen_t pos = 0;
  for (R_xlen_t i = 0; i < nr; ++i) {
    const double ri = r[i];
    double x = x0;
    for (int k = 0; k < n_transient; ++k) x = ri * x * (1.0 - x);
    for (int k = 0; k < n_keep; ++k) {
      x = ri * x * (1.0 - x);
      rr[pos] = ri;
      xx[pos] = x;
      ++pos;
    }
  }
  return List::create(_["r"] = rr, _["x"] = xx);
}

// Lyapunov exponent of the logistic map for each r:
//   lambda(r) = lim (1/n) sum log |f'(x_k)|,  f'(x) = r (1 - 2x)
// [[Rcpp::export]]
NumericVector logistic_lyapunov_cpp(NumericVector r, double x0, int n_transient,
                                    int n) {
  const R_xlen_t nr = r.size();
  NumericVector out(nr);
  for (R_xlen_t i = 0; i < nr; ++i) {
    const double ri = r[i];
    double x = x0;
    for (int k = 0; k < n_transient; ++k) x = ri * x * (1.0 - x);
    double acc = 0.0;
    for (int k = 0; k < n; ++k) {
      double d = std::fabs(ri * (1.0 - 2.0 * x));
      // A superstable orbit hits x = 1/2 exactly and gives log(0); clamp so
      // the estimate stays finite (and very negative).
      if (d < 1e-300) d = 1e-300;
      acc += std::log(d);
      x = ri * x * (1.0 - x);
    }
    out[i] = acc / static_cast<double>(n);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Henon map  x_{n+1} = 1 - a x_n^2 + y_n,  y_{n+1} = b x_n
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix henon_cpp(double x0, double y0, double a, double b, int n,
                        int n_transient) {
  if (n < 1) stop("n must be at least 1");
  double x = x0, y = y0;
  NumericMatrix out(n, 2);
  for (int k = 0; k < n_transient; ++k) {
    const double xn = 1.0 - a * x * x + y;
    y = b * x;
    x = xn;
    if (!std::isfinite(x) || std::fabs(x) > 1e6) {
      std::fill(out.begin(), out.end(), NA_REAL);
      return out;
    }
  }
  for (int k = 0; k < n; ++k) {
    out(k, 0) = x;
    out(k, 1) = y;
    const double xn = 1.0 - a * x * x + y;
    y = b * x;
    x = xn;
    if (!std::isfinite(x) || std::fabs(x) > 1e6) {
      // Orbit escaped to infinity: fill the remainder with NA and stop.
      for (int j = k + 1; j < n; ++j) {
        out(j, 0) = NA_REAL;
        out(j, 1) = NA_REAL;
      }
      break;
    }
  }
  return out;
}

// [[Rcpp::export]]
List henon_bifurcation_cpp(NumericVector a, double b, double x0, double y0,
                           int n_transient, int n_keep) {
  const R_xlen_t na = a.size();
  const R_xlen_t total = na * static_cast<R_xlen_t>(n_keep);
  NumericVector aa(total), xx(total);
  R_xlen_t pos = 0;
  for (R_xlen_t i = 0; i < na; ++i) {
    const double ai = a[i];
    double x = x0, y = y0;
    bool escaped = false;
    for (int k = 0; k < n_transient; ++k) {
      const double xn = 1.0 - ai * x * x + y;
      y = b * x;
      x = xn;
      if (!std::isfinite(x) || std::fabs(x) > 1e6) { escaped = true; break; }
    }
    for (int k = 0; k < n_keep; ++k) {
      aa[pos] = ai;
      if (escaped) {
        xx[pos] = NA_REAL;
      } else {
        const double xn = 1.0 - ai * x * x + y;
        y = b * x;
        x = xn;
        if (!std::isfinite(x) || std::fabs(x) > 1e6) { escaped = true; xx[pos] = NA_REAL; }
        else xx[pos] = x;
      }
      ++pos;
    }
  }
  return List::create(_["a"] = aa, _["x"] = xx);
}

// ---------------------------------------------------------------------------
// Ikeda map (laser cavity model), u = 0.918 gives the classic attractor
//   t = 0.4 - 6 / (1 + x^2 + y^2)
//   x_{n+1} = 1 + u (x cos t - y sin t)
//   y_{n+1} =     u (x sin t + y cos t)
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix ikeda_cpp(double x0, double y0, double u, int n, int n_transient) {
  if (n < 1) stop("n must be at least 1");
  double x = x0, y = y0;
  auto step = [u](double& x, double& y) {
    const double t = 0.4 - 6.0 / (1.0 + x * x + y * y);
    const double ct = std::cos(t), st = std::sin(t);
    const double xn = 1.0 + u * (x * ct - y * st);
    const double yn = u * (x * st + y * ct);
    x = xn;
    y = yn;
  };
  for (int k = 0; k < n_transient; ++k) step(x, y);
  NumericMatrix out(n, 2);
  for (int k = 0; k < n; ++k) {
    out(k, 0) = x;
    out(k, 1) = y;
    step(x, y);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Chirikov standard map on the torus [0, 2pi)^2
//   p_{n+1}     = p_n + K sin(theta_n)
//   theta_{n+1} = theta_n + p_{n+1}
// Several initial conditions are iterated so islands and chaotic seas can
// be seen side by side. Output columns: orbit id, theta, p.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix standard_map_cpp(NumericVector theta0, NumericVector p0, double K,
                               int n) {
  if (theta0.size() != p0.size()) stop("theta0 and p0 must have the same length");
  if (n < 1) stop("n must be at least 1");
  const R_xlen_t m = theta0.size();
  NumericMatrix out(m * static_cast<R_xlen_t>(n), 3);
  R_xlen_t pos = 0;
  for (R_xlen_t j = 0; j < m; ++j) {
    double th = wrap_2pi(theta0[j]);
    double p = wrap_2pi(p0[j]);
    for (int k = 0; k < n; ++k) {
      out(pos, 0) = static_cast<double>(j + 1);
      out(pos, 1) = th;
      out(pos, 2) = p;
      ++pos;
      p = wrap_2pi(p + K * std::sin(th));
      th = wrap_2pi(th + p);
    }
  }
  return out;
}
