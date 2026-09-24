// Recurrence plots and Recurrence Quantification Analysis (RQA).
//
// A recurrence matrix is R(i, j) = 1 when ||x_i - x_j|| <= eps for state
// vectors x_i of an (optionally delay-embedded) time series. RQA measures
// summarise the diagonal and vertical line structures of that matrix
// (Marwan et al., 2007, Physics Reports 438, 237-329). The definitions below
// follow that paper and the CRP Toolbox conventions:
//
//   RR   recurrence rate, density of recurrence points
//   DET  determinism, fraction of points on diagonal lines of length >= lmin
//   L    average diagonal line length
//   Lmax longest diagonal line, DIV = 1 / Lmax
//   ENTR Shannon entropy of the diagonal line length distribution
//   LAM  laminarity, fraction of points on vertical lines of length >= vmin
//   TT   trapping time, average vertical line length
//   Vmax longest vertical line
//
// The Theiler window removes cells with |i - j| < theiler from every
// computation (theiler = 1 drops the main diagonal, the line of identity).

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

// Row-major copy so that a state vector is contiguous in memory.
std::vector<double> to_row_major(const NumericMatrix& X) {
  const int N = X.nrow(), m = X.ncol();
  std::vector<double> out(static_cast<size_t>(N) * m);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < m; ++j) out[static_cast<size_t>(i) * m + j] = X(i, j);
  return out;
}

// norm: 1 = Euclidean (L2), 2 = maximum (L-infinity), 3 = Manhattan (L1)
inline double dist_pair(const double* a, const double* b, int m, int norm) {
  if (norm == 2) {
    double d = 0.0;
    for (int k = 0; k < m; ++k) d = std::max(d, std::fabs(a[k] - b[k]));
    return d;
  }
  if (norm == 3) {
    double d = 0.0;
    for (int k = 0; k < m; ++k) d += std::fabs(a[k] - b[k]);
    return d;
  }
  double s = 0.0;
  for (int k = 0; k < m; ++k) {
    const double diff = a[k] - b[k];
    s += diff * diff;
  }
  return std::sqrt(s);
}

// Scan a run-length histogram along one line of the matrix.
// (i, j) is the start cell, (di, dj) the step.
inline void scan_line(const std::vector<unsigned char>& R, size_t N, long i,
                      long j, long di, long dj, std::vector<double>& hist) {
  long len = 0;
  const long n = static_cast<long>(N);
  while (i >= 0 && j >= 0 && i < n && j < n) {
    if (R[static_cast<size_t>(i) * N + static_cast<size_t>(j)]) {
      ++len;
    } else if (len > 0) {
      hist[static_cast<size_t>(len)] += 1.0;
      len = 0;
    }
    i += di;
    j += dj;
  }
  if (len > 0) hist[static_cast<size_t>(len)] += 1.0;
}

}  // namespace

// [[Rcpp::export]]
NumericMatrix embed_cpp(NumericVector x, int m, int tau) {
  const R_xlen_t N = x.size();
  if (m < 1) stop("m must be at least 1");
  if (tau < 1) stop("tau must be at least 1");
  const R_xlen_t n_rows = N - static_cast<R_xlen_t>(m - 1) * tau;
  if (n_rows < 1) stop("series too short for this embedding (m = %d, tau = %d)", m, tau);
  NumericMatrix out(n_rows, m);
  for (int j = 0; j < m; ++j) {
    const R_xlen_t off = static_cast<R_xlen_t>(j) * tau;
    for (R_xlen_t i = 0; i < n_rows; ++i) out(i, j) = x[i + off];
  }
  return out;
}

// [[Rcpp::export]]
LogicalMatrix recurrence_matrix_cpp(NumericMatrix X, double eps, int norm) {
  const int N = X.nrow(), m = X.ncol();
  const std::vector<double> A = to_row_major(X);
  LogicalMatrix R(N, N);
  for (int i = 0; i < N; ++i) {
    R(i, i) = TRUE;
    const double* a = &A[static_cast<size_t>(i) * m];
    for (int j = i + 1; j < N; ++j) {
      const bool rec = dist_pair(a, &A[static_cast<size_t>(j) * m], m, norm) <= eps;
      R(i, j) = rec;
      R(j, i) = rec;
    }
  }
  return R;
}

// Sample (or enumerate, when small) pairwise distances. Used on the R side
// to pick eps for a target recurrence rate.
// [[Rcpp::export]]
NumericVector sample_distances_cpp(NumericMatrix X, int n_pairs, int norm) {
  const int N = X.nrow(), m = X.ncol();
  if (N < 2) stop("need at least two state vectors");
  const std::vector<double> A = to_row_major(X);
  const double n_all = static_cast<double>(N) * (N - 1) / 2.0;

  if (n_all <= static_cast<double>(n_pairs)) {
    NumericVector out(static_cast<R_xlen_t>(n_all));
    R_xlen_t pos = 0;
    for (int i = 0; i < N; ++i)
      for (int j = i + 1; j < N; ++j)
        out[pos++] = dist_pair(&A[static_cast<size_t>(i) * m],
                               &A[static_cast<size_t>(j) * m], m, norm);
    return out;
  }

  RNGScope scope;  // uses R's RNG, so set.seed() makes this reproducible
  NumericVector out(n_pairs);
  for (int k = 0; k < n_pairs; ++k) {
    int i = static_cast<int>(unif_rand() * N);
    int j = static_cast<int>(unif_rand() * N);
    while (j == i) j = static_cast<int>(unif_rand() * N);
    if (i >= N) i = N - 1;
    if (j >= N) j = N - 1;
    out[k] = dist_pair(&A[static_cast<size_t>(i) * m],
                       &A[static_cast<size_t>(j) * m], m, norm);
  }
  return out;
}

// [[Rcpp::export]]
List rqa_cpp(NumericMatrix X, double eps, int norm, int lmin, int vmin, int theiler) {
  const int N = X.nrow(), m = X.ncol();
  if (N < 2) stop("need at least two state vectors");
  if (lmin < 1 || vmin < 1) stop("lmin and vmin must be at least 1");
  if (theiler < 0) stop("theiler must be non-negative");

  const size_t NN = static_cast<size_t>(N);
  const std::vector<double> A = to_row_major(X);

  // Build the recurrence matrix with the Theiler window already applied.
  std::vector<unsigned char> R(NN * NN, 0);
  double n_rec = 0.0;
  if (theiler == 0) {
    for (size_t i = 0; i < NN; ++i) R[i * NN + i] = 1;
    n_rec += static_cast<double>(N);
  }
  const int start_off = std::max(theiler, 1);
  for (int i = 0; i < N; ++i) {
    const double* a = &A[static_cast<size_t>(i) * m];
    for (int j = i + start_off; j < N; ++j) {
      if (dist_pair(a, &A[static_cast<size_t>(j) * m], m, norm) <= eps) {
        R[static_cast<size_t>(i) * NN + j] = 1;
        R[static_cast<size_t>(j) * NN + i] = 1;
        n_rec += 2.0;
      }
    }
  }

  // Cells excluded by the Theiler window: |i - j| < theiler.
  double excluded = 0.0;
  if (theiler >= 1) excluded += static_cast<double>(N);
  for (int d = 1; d < theiler; ++d) excluded += 2.0 * (N - d);
  const double admissible = static_cast<double>(N) * N - excluded;

  // Diagonal line histogram (both triangles, |offset| >= theiler).
  std::vector<double> diag_hist(NN + 1, 0.0);
  for (long d = -(N - 1); d <= N - 1; ++d) {
    if (std::labs(d) < theiler) continue;
    const long i0 = d < 0 ? -d : 0;
    const long j0 = d > 0 ? d : 0;
    scan_line(R, NN, i0, j0, 1, 1, diag_hist);
  }

  // Vertical line histogram (columns; the matrix is symmetric so this equals
  // the horizontal one).
  std::vector<double> vert_hist(NN + 1, 0.0);
  for (long j = 0; j < N; ++j) scan_line(R, NN, 0, j, 1, 0, vert_hist);

  // Measures from the diagonal structures.
  double det_num = 0.0, l_count = 0.0, lmax = 0.0;
  for (size_t l = static_cast<size_t>(lmin); l <= NN; ++l) {
    const double P = diag_hist[l];
    if (P > 0.0) {
      det_num += static_cast<double>(l) * P;
      l_count += P;
      lmax = static_cast<double>(l);
    }
  }
  double entr = 0.0;
  if (l_count > 0.0) {
    for (size_t l = static_cast<size_t>(lmin); l <= NN; ++l) {
      const double P = diag_hist[l];
      if (P > 0.0) {
        const double p = P / l_count;
        entr -= p * std::log(p);
      }
    }
  }

  // Measures from the vertical structures.
  double lam_num = 0.0, v_count = 0.0, vmax = 0.0;
  for (size_t v = static_cast<size_t>(vmin); v <= NN; ++v) {
    const double P = vert_hist[v];
    if (P > 0.0) {
      lam_num += static_cast<double>(v) * P;
      v_count += P;
      vmax = static_cast<double>(v);
    }
  }

  const double RR = admissible > 0.0 ? n_rec / admissible : NA_REAL;
  const double DET = n_rec > 0.0 ? det_num / n_rec : NA_REAL;
  const double L = l_count > 0.0 ? det_num / l_count : NA_REAL;
  const double LAM = n_rec > 0.0 ? lam_num / n_rec : NA_REAL;
  const double TT = v_count > 0.0 ? lam_num / v_count : NA_REAL;
  const double DIV = lmax > 0.0 ? 1.0 / lmax : NA_REAL;
  const double RATIO = (RR > 0.0 && n_rec > 0.0) ? DET / RR : NA_REAL;

  // Drop the unused length-0 slot so hist[l] is the count of lines of length l.
  NumericVector dh(diag_hist.begin() + 1, diag_hist.end());
  NumericVector vh(vert_hist.begin() + 1, vert_hist.end());

  return List::create(
      _["RR"] = RR, _["DET"] = DET, _["L"] = L, _["Lmax"] = lmax, _["DIV"] = DIV,
      _["ENTR"] = entr, _["LAM"] = LAM, _["TT"] = TT, _["Vmax"] = vmax,
      _["RATIO"] = RATIO, _["n_recurrences"] = n_rec, _["n_states"] = N,
      _["diag_hist"] = dh, _["vert_hist"] = vh);
}
