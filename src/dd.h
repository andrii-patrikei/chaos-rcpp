// Double-double arithmetic: a number is the unevaluated sum hi + lo of two
// IEEE doubles, giving about 106 bits (32 decimal digits) of precision with
// nothing but ordinary double operations. Built from the classical
// error-free transformations (Knuth's TwoSum, Dekker's TwoProd via fma),
// following Hida, Li and Bailey (2001), "Algorithms for quad-double
// precision floating point arithmetic".
//
// The type provides exactly the operators the templated RK4 kernel needs, so
// the same integrator code that runs in float and double also runs here.
// It is much slower than double (roughly 50-90x in this integrator, because
// every double-double addition is a chain of about twenty dependent
// floating-point operations and the code becomes latency bound) but still
// about 2000x faster than a pure R loop over Rmpfr numbers, which makes it a
// practical high-precision reference for long runs. The Rmpfr path on the R
// side remains the gold standard for validating it.

#ifndef CHAOSRCPP_DD_H
#define CHAOSRCPP_DD_H

#include <cmath>

namespace chaos {

struct dd {
  double hi, lo;

  dd() : hi(0.0), lo(0.0) {}
  dd(double h) : hi(h), lo(0.0) {}  // NOLINT: implicit on purpose, mirrors double
  dd(double h, double l) : hi(h), lo(l) {}

  explicit operator double() const { return hi + lo; }
};

// s + e == a + b exactly
inline dd two_sum(double a, double b) {
  const double s = a + b;
  const double bb = s - a;
  const double e = (a - (s - bb)) + (b - bb);
  return dd(s, e);
}

// requires |a| >= |b|
inline dd quick_two_sum(double a, double b) {
  const double s = a + b;
  const double e = b - (s - a);
  return dd(s, e);
}

// p + e == a * b exactly
inline dd two_prod(double a, double b) {
  const double p = a * b;
  const double e = std::fma(a, b, -p);
  return dd(p, e);
}

inline dd operator+(const dd& a, const dd& b) {
  dd s = two_sum(a.hi, b.hi);
  const dd t = two_sum(a.lo, b.lo);
  s.lo += t.hi;
  s = quick_two_sum(s.hi, s.lo);
  s.lo += t.lo;
  return quick_two_sum(s.hi, s.lo);
}

inline dd operator-(const dd& a) { return dd(-a.hi, -a.lo); }

inline dd operator-(const dd& a, const dd& b) { return a + (-b); }

inline dd operator*(const dd& a, const dd& b) {
  dd p = two_prod(a.hi, b.hi);
  p.lo += a.hi * b.lo + a.lo * b.hi;
  return quick_two_sum(p.hi, p.lo);
}

inline dd operator/(const dd& a, const dd& b) {
  // Long division with three correction steps (Hida et al., accurate div).
  const double q1 = a.hi / b.hi;
  dd r = a - b * dd(q1);
  const double q2 = r.hi / b.hi;
  r = r - b * dd(q2);
  const double q3 = r.hi / b.hi;
  const dd q = quick_two_sum(q1, q2);
  return q + dd(q3);
}

inline dd& operator+=(dd& a, const dd& b) { a = a + b; return a; }
inline dd& operator-=(dd& a, const dd& b) { a = a - b; return a; }
inline dd& operator*=(dd& a, const dd& b) { a = a * b; return a; }

}  // namespace chaos

#endif  // CHAOSRCPP_DD_H
