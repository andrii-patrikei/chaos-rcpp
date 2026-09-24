// Continuous-time chaotic flows integrated with a classical 4th-order
// Runge-Kutta scheme. Everything is templated on the floating-point type T,
// so the very same code path runs in single precision (float, fp32), double
// precision (double, fp64) and double-double (chaos::dd, about 106 bits).
// The R side adds an Rmpfr run that uses the same RK4 formulas in arbitrary
// precision as the gold-standard reference.
//
// precision codes used by the exported functions: 0 = fp64, 1 = fp32, 2 = dd

#include <Rcpp.h>
#include <cmath>
#include <vector>

#include "dd.h"

using namespace Rcpp;
using chaos::dd;

// ---------------------------------------------------------------------------
// Right-hand sides
// ---------------------------------------------------------------------------

template <typename T>
struct Lorenz {
  T sigma, rho, beta;
  static constexpr int dim = 3;

  inline void operator()(const T* x, T* dx) const {
    dx[0] = sigma * (x[1] - x[0]);
    dx[1] = x[0] * (rho - x[2]) - x[1];
    dx[2] = x[0] * x[1] - beta * x[2];
  }

  // Jacobian-vector product J(x) v, used for the tangent dynamics.
  inline void jac(const T* x, const T* v, T* dv) const {
    dv[0] = sigma * (v[1] - v[0]);
    dv[1] = (rho - x[2]) * v[0] - v[1] - x[0] * v[2];
    dv[2] = x[1] * v[0] + x[0] * v[1] - beta * v[2];
  }
};

template <typename T>
struct Roessler {
  T a, b, c;
  static constexpr int dim = 3;

  inline void operator()(const T* x, T* dx) const {
    dx[0] = -x[1] - x[2];
    dx[1] = x[0] + a * x[1];
    dx[2] = b + x[2] * (x[0] - c);
  }

  inline void jac(const T* x, const T* v, T* dv) const {
    dv[0] = -v[1] - v[2];
    dv[1] = v[0] + a * v[1];
    dv[2] = x[2] * v[0] + (x[0] - c) * v[2];
  }
};

// State + D tangent vectors, integrated as one larger system.
template <typename T, typename Sys>
struct TangentSystem {
  const Sys& f;
  static constexpr int D = Sys::dim;
  static constexpr int dim = D + D * D;

  explicit TangentSystem(const Sys& sys) : f(sys) {}

  inline void operator()(const T* s, T* ds) const {
    f(s, ds);
    for (int j = 0; j < D; ++j) {
      f.jac(s, s + D + j * D, ds + D + j * D);
    }
  }
};

// ---------------------------------------------------------------------------
// Classical RK4 step, generic in precision T and dimension D
// ---------------------------------------------------------------------------

template <typename T, int D, typename Sys>
inline void rk4_step(const Sys& f, T* x, const T h) {
  T k1[D], k2[D], k3[D], k4[D], tmp[D];
  const T half = T(0.5);
  const T two = T(2);
  const T sixth = T(1) / T(6);

  f(x, k1);
  for (int i = 0; i < D; ++i) tmp[i] = x[i] + half * h * k1[i];
  f(tmp, k2);
  for (int i = 0; i < D; ++i) tmp[i] = x[i] + half * h * k2[i];
  f(tmp, k3);
  for (int i = 0; i < D; ++i) tmp[i] = x[i] + h * k3[i];
  f(tmp, k4);
  for (int i = 0; i < D; ++i) {
    x[i] += sixth * h * (k1[i] + two * k2[i] + two * k3[i] + k4[i]);
  }
}

// How a state of type T is written to the output matrix. float and double
// occupy D columns; a double-double state occupies 2D columns (all hi parts,
// then all lo parts) so that no precision is lost on the way to R.
template <typename T>
struct out_traits {
  static int cols(int D) { return D; }
  static void store(NumericMatrix& out, int row, int i, int, const T& v) {
    out(row, i) = static_cast<double>(v);
  }
};

template <>
struct out_traits<dd> {
  static int cols(int D) { return 2 * D; }
  static void store(NumericMatrix& out, int row, int i, int D, const dd& v) {
    out(row, i) = v.hi;
    out(row, D + i) = v.lo;
  }
};

// Integrate n steps of size dt, storing every `thin`-th state.
template <typename T, typename Sys>
NumericMatrix integrate_flow(const Sys& f, const NumericVector& x0,
                             const double dt, const int n, const int thin) {
  constexpr int D = Sys::dim;
  if (x0.size() != D) stop("x0 must have length %d", D);
  if (n < 1) stop("n must be at least 1");
  if (thin < 1) stop("thin must be at least 1");

  const int n_out = n / thin + 1;
  NumericMatrix out(n_out, out_traits<T>::cols(D));

  T x[D];
  for (int i = 0; i < D; ++i) x[i] = static_cast<T>(x0[i]);
  const T h = static_cast<T>(dt);

  for (int i = 0; i < D; ++i) out_traits<T>::store(out, 0, i, D, x[i]);

  int row = 1;
  for (int k = 1; k <= n; ++k) {
    rk4_step<T, D>(f, x, h);
    if (k % thin == 0 && row < n_out) {
      for (int i = 0; i < D; ++i) out_traits<T>::store(out, row, i, D, x[i]);
      ++row;
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// Lyapunov spectrum by Gram-Schmidt reorthonormalisation (Benettin et al.)
// ---------------------------------------------------------------------------

template <typename Sys>
NumericVector lyapunov_spectrum_flow(const Sys& f, const NumericVector& x0,
                                     const double dt, const int n,
                                     const int n_transient) {
  constexpr int D = Sys::dim;
  using TS = TangentSystem<double, Sys>;
  constexpr int DD = TS::dim;
  if (x0.size() != D) stop("x0 must have length %d", D);

  TS ts(f);
  double s[DD];
  for (int i = 0; i < D; ++i) s[i] = x0[i];

  // Let the state settle onto the attractor first.
  for (int k = 0; k < n_transient; ++k) rk4_step<double, D>(f, s, dt);

  // Tangent vectors start as the identity (column j lives at s + D + j*D).
  for (int i = 0; i < D * D; ++i) s[D + i] = 0.0;
  for (int j = 0; j < D; ++j) s[D + j * D + j] = 1.0;

  std::vector<double> sum(D, 0.0);

  for (int k = 0; k < n; ++k) {
    rk4_step<double, DD>(ts, s, dt);

    // Modified Gram-Schmidt on the D tangent vectors.
    for (int j = 0; j < D; ++j) {
      double* vj = s + D + j * D;
      for (int p = 0; p < j; ++p) {
        const double* vp = s + D + p * D;
        double dot = 0.0;
        for (int i = 0; i < D; ++i) dot += vj[i] * vp[i];
        for (int i = 0; i < D; ++i) vj[i] -= dot * vp[i];
      }
      double norm = 0.0;
      for (int i = 0; i < D; ++i) norm += vj[i] * vj[i];
      norm = std::sqrt(norm);
      sum[j] += std::log(norm);
      for (int i = 0; i < D; ++i) vj[i] /= norm;
    }
  }

  NumericVector out(D);
  const double total_time = static_cast<double>(n) * dt;
  for (int j = 0; j < D; ++j) out[j] = sum[j] / total_time;
  return out;
}

// ---------------------------------------------------------------------------
// Exported entry points
// ---------------------------------------------------------------------------

template <template <typename> class Sys, typename... Args>
NumericMatrix dispatch_flow(int precision, const NumericVector& x0, double dt,
                            int n, int thin, Args... p) {
  switch (precision) {
    case 0: {
      Sys<double> f{static_cast<double>(p)...};
      return integrate_flow<double>(f, x0, dt, n, thin);
    }
    case 1: {
      Sys<float> f{static_cast<float>(p)...};
      return integrate_flow<float>(f, x0, dt, n, thin);
    }
    case 2: {
      Sys<dd> f{dd(p)...};
      return integrate_flow<dd>(f, x0, dt, n, thin);
    }
    default:
      stop("unknown precision code %d", precision);
  }
}

// [[Rcpp::export]]
NumericMatrix lorenz_cpp(NumericVector x0, double sigma, double rho, double beta,
                         double dt, int n, int thin, int precision) {
  return dispatch_flow<Lorenz>(precision, x0, dt, n, thin, sigma, rho, beta);
}

// [[Rcpp::export]]
NumericMatrix roessler_cpp(NumericVector x0, double a, double b, double c,
                           double dt, int n, int thin, int precision) {
  return dispatch_flow<Roessler>(precision, x0, dt, n, thin, a, b, c);
}

// [[Rcpp::export]]
NumericVector lorenz_lyapunov_cpp(NumericVector x0, double sigma, double rho,
                                  double beta, double dt, int n, int n_transient) {
  Lorenz<double> f{sigma, rho, beta};
  return lyapunov_spectrum_flow(f, x0, dt, n, n_transient);
}

// [[Rcpp::export]]
NumericVector roessler_lyapunov_cpp(NumericVector x0, double a, double b, double c,
                                    double dt, int n, int n_transient) {
  Roessler<double> f{a, b, c};
  return lyapunov_spectrum_flow(f, x0, dt, n, n_transient);
}

// Tiny self-test of the double-double arithmetic: returns the hi and lo parts
// of a few operations so the R side can compare them with Rmpfr.
// [[Rcpp::export]]
NumericVector dd_selftest_cpp(double a, double b) {
  const dd x(a), y(b);
  const dd s = x + y, d = x - y, p = x * y, q = x / y;
  return NumericVector::create(s.hi, s.lo, d.hi, d.lo, p.hi, p.lo, q.hi, q.lo);
}
