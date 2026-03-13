#include "XOMultinom.h"

// [[Rcpp::depends("RcppArmadillo")]]

// Note: RcppExport is an alias for extern "C"

static void vstep(std::vector<double>& v, std::vector<double>& w,
                  int n, double p, int b, int a) {
  const double q  = 1.0 - p;
  const double pq = (q > 1e-300) ? p / q : 0.0;

  std::fill(w.begin(), w.end(), 0.0);

  // Edge case: p = 1  =>  Binom(n-i,1) concentrates all mass at d = n-i
  if (p > 1.0 - 1e-14) {
    for (int i = 0; i <= n; ++i) {
      if (v[i] == 0.0) continue;
      const int d = n - i;
      if (d >= b && d <= a) w[n] += v[i];
    }
    std::swap(v, w);
    return;
  }

  // Precompute q^(n-i) for i = 0 .. n  in O(n)
  std::vector<double> qpow(n + 1);
  qpow[n] = 1.0;
  for (int i = n - 1; i >= 0; --i) qpow[i] = qpow[i + 1] * q;

  for (int i = 0; i <= n; ++i) {
    if (v[i] == 0.0) continue;
    const int max_d = n - i;
    if (b > max_d) continue;          // can't place >= b balls

    // Start at d=0: f = Binom(max_d, p, 0) = q^max_d
    double f = qpow[i];

    // Advance recurrence from d=0 up to d=b
    //   Binom(max_d, p, d) = Binom(max_d, p, d-1) * (p/q) * (max_d-d+1)/d
    for (int d = 1; d <= b; ++d)
      f *= pq * (max_d - d + 1.0) / d;
    // f = Binom(max_d, p, b)

    const int  d_hi = std::min(a, max_d);
    const double vi = v[i];
    for (int d = b; d <= d_hi; ++d) {
      w[i + d] += vi * f;
      if (d < d_hi)
        f *= pq * (max_d - d) / (d + 1.0);
    }
  }
  std::swap(v, w);
}

static inline bool max_convert(double c, int n, int& ci) {
  if (std::isnan(c))            { ci = -1; return false; }
  if (c < 0.0)                  { ci = -1; return true;  }
  if (c >= static_cast<double>(n)) { ci = n; return true;  }
  ci = static_cast<int>(std::floor(c));
  return true;
}

static inline bool min_convert(double c, int n, int& ci) {
  if (std::isnan(c))               { ci =  0; return false; }
  if (c <= 0.0)                    { ci =  0; return true;  }
  if (c > static_cast<double>(n)) { ci = n + 1; return true; }
  ci = static_cast<int>(std::ceil(c));
  return true;
}

// [[Rcpp::export]]
double prob_max_leq(int n, const std::vector<double>& pi, double c) {
  int ci;
  if (!max_convert(c, n, ci)) return 0.0;   // NaN → 0
  if (ci >= n) return 1.0;
  if (ci <  0) return 0.0;

  const int m = static_cast<int>(pi.size());

  // suffix_sum[k] = sum_{j=k}^{m-1} pi[j]
  std::vector<double> suf(m + 1, 0.0);
  for (int k = m - 1; k >= 0; --k) suf[k] = suf[k + 1] + pi[k];

  std::vector<double> v(n + 1, 0.0), w(n + 1);

  // Initiating row vector: culled first row of Q_1  (pi_0* = pi[0])
  {
    const double p  = pi[0];
    const double q  = 1.0 - p;
    const double pq = (q > 1e-300) ? p / q : 0.0;
    if (p > 1.0 - 1e-14) {
      if (n <= ci) v[n] = 1.0;
    } else {
      double f = std::pow(q, static_cast<double>(n));   // Binom(n,p,0)
      v[0] = f;
      for (int j = 1; j <= std::min(ci, n); ++j) {
        f *= pq * (n - j + 1.0) / j;
        v[j] = f;
      }
      // v[j] = 0 for j > ci  (initialized to 0)
    }
  }

  // Propagate through urns k = 1 .. m-2  (0-indexed)
  for (int k = 1; k <= m - 2; ++k)
    vstep(v, w, n, pi[k] / suf[k], 0, ci);

  // Final summation: last urn gets n_m = n - j balls.
  // Require n_m <= ci  =>  j >= n - ci
  double result = 0.0;
  for (int j = std::max(0, n - ci); j <= n; ++j) result += v[j];
  return std::min(1.0, result);
}

// [[Rcpp::export]]
double prob_min_geq(int n, const std::vector<double>& pi, double c) {
  int ci;
  if (!min_convert(c, n, ci)) return 1.0;   // NaN --> 1  (conservative)
  if (ci <= 0)                return 1.0;
  if (static_cast<long long>(ci) * pi.size() > static_cast<long long>(n)) return 0.0;

  const int m = static_cast<int>(pi.size());

  std::vector<double> suf(m + 1, 0.0);
  for (int k = m - 1; k >= 0; --k) suf[k] = suf[k + 1] + pi[k];

  std::vector<double> v(n + 1, 0.0), w(n + 1);

  // Initiating row vector: culled first row of Q_1  (keep d >= ci)
  {
    const double p  = pi[0];
    const double q  = 1.0 - p;
    const double pq = (q > 1e-300) ? p / q : 0.0;
    if (p > 1.0 - 1e-14) {
      if (n >= ci) v[n] = 1.0;
    } else {
      double f = std::pow(q, static_cast<double>(n));   // Binom(n,p,0)
      // Advance to j = ci
      for (int j = 1; j <= std::min(ci, n); ++j)
        f *= pq * (n - j + 1.0) / j;
      // f = Binom(n, p, ci)
      for (int j = ci; j <= n; ++j) {
        v[j] = f;
        if (j < n) f *= pq * (n - j) / (j + 1.0);
      }
    }
  }

  // Propagate through urns k = 1 .. m-2
  for (int k = 1; k <= m - 2; ++k)
    vstep(v, w, n, pi[k] / suf[k], ci, n);

  // Final summation: require n_m = n - j >= ci  =>  j <= n - ci
  double result = 0.0;
  for (int j = 0; j <= n - ci; ++j) result += v[j];
  return std::min(1.0, result);
}

// Convenience wrappers
// [[Rcpp::export]]
double prob_max_gt(int n, const std::vector<double>& pi, double c) {
  return std::max(0.0, 1.0 - prob_max_leq(n, pi, c));
}

// [[Rcpp::export]]
double prob_min_lt(int n, const std::vector<double>& pi, double c) {
  return std::max(0.0, 1.0 - prob_min_geq(n, pi, c));
}

// [[Rcpp::export]]
double prob_min_leq(int n, const std::vector<double>& pi, double c) {
  if (std::isnan(c)) return 0.0;
  if (c < 0.0)       return 0.0;   // min is always >= 0
  if (c >= static_cast<double>(n)) return 1.0;
  // floor(c) is now safe to compute
  return std::max(0.0, 1.0 - prob_min_geq(n, pi, std::floor(c) + 1.0));
}

// [[Rcpp::export]]
double prob_joint(int n, const std::vector<double>& pi, double a, double b) {
  // Convert bounds
  int ai, bi;
  if (!max_convert(a, n, ai)) return 0.0;
  if (!min_convert(b, n, bi)) return 0.0;

  if (ai >= n && bi <= 0) return 1.0;
  if (ai < bi || ai < 0)  return 0.0;
  if (static_cast<long long>(bi) * pi.size() > static_cast<long long>(n)) return 0.0;

  const int m = static_cast<int>(pi.size());

  std::vector<double> suf(m + 1, 0.0);
  for (int k = m - 1; k >= 0; --k) suf[k] = suf[k + 1] + pi[k];

  std::vector<double> v(n + 1, 0.0), w(n + 1);

  // Initiating vector: keep bi <= d <= ai
  {
    const double p  = pi[0];
    const double q  = 1.0 - p;
    const double pq = (q > 1e-300) ? p / q : 0.0;
    const int lo = std::max(0, bi), hi = std::min(ai, n);
    if (p > 1.0 - 1e-14) {
      if (n >= lo && n <= hi) v[n] = 1.0;
    } else {
      double f = std::pow(q, static_cast<double>(n));
      for (int j = 1; j <= lo; ++j)
        f *= pq * (n - j + 1.0) / j;
      for (int j = lo; j <= hi; ++j) {
        v[j] = f;
        if (j < hi) f *= pq * (n - j) / (j + 1.0);
      }
    }
  }

  for (int k = 1; k <= m - 2; ++k)
    vstep(v, w, n, pi[k] / suf[k], bi, ai);

  // Final sum: bi <= n - j <= ai  =>  n-ai <= j <= n-bi
  double result = 0.0;
  for (int j = std::max(0, n - ai); j <= std::min(n, n - bi); ++j)
    result += v[j];
  return std::min(1.0, result);
}

// PMF: P(max n_k = c)
// [[Rcpp::export]]
double prob_max_eq(int n, const std::vector<double>& pi, double c) {
  // Non-integer, NaN, or Inf → PMF is 0
  if (!std::isfinite(c) || c != std::floor(c)) return 0.0;
  const int ci = static_cast<int>(c);
  if (ci < 0 || ci > n) return 0.0;
  double hi = prob_max_leq(n, pi, static_cast<double>(ci));
  double lo = (ci > 0) ? prob_max_leq(n, pi, static_cast<double>(ci - 1)) : 0.0;
  return std::max(0.0, hi - lo);
}

// PMF: P(min n_k = c)
// [[Rcpp::export]]
double prob_min_eq(int n, const std::vector<double>& pi, double c) {
  if (!std::isfinite(c) || c != std::floor(c)) return 0.0;
  const int ci = static_cast<int>(c);
  if (ci < 0 || ci > n) return 0.0;
  double hi = prob_min_geq(n, pi, static_cast<double>(ci));
  double lo = prob_min_geq(n, pi, static_cast<double>(ci + 1));
  return std::max(0.0, hi - lo);
}

// PMF vector for max over c in [c_lo, c_hi]
// Chains adjacent CDF calls: each curr becomes the next prev, so the
// boundary between successive PMF values is computed only once.
// [[Rcpp::export]]
Rcpp::NumericVector pmf_max_range(int n, const std::vector<double>& pi,
                                  double c_lo, double c_hi) {
  if (!std::isfinite(c_lo) || !std::isfinite(c_hi))
    return Rcpp::NumericVector(0);

  // First integer >= c_lo, last integer <= c_hi
  const int lo = std::max(0, static_cast<int>(std::ceil(c_lo)));
  const int hi = std::min(n, static_cast<int>(std::floor(c_hi)));
  if (lo > hi) return Rcpp::NumericVector(0);

  const int len = hi - lo + 1;
  Rcpp::NumericVector out(len);
  double prev = (lo > 0) ? prob_max_leq(n, pi, static_cast<double>(lo - 1)) : 0.0;
  for (int i = 0; i < len; ++i) {
    double curr = prob_max_leq(n, pi, static_cast<double>(lo + i));
    out[i] = std::max(0.0, curr - prev);
    prev = curr;
  }
  return out;
}

// PMF vector for min over c in [c_lo, c_hi]
// [[Rcpp::export]]
Rcpp::NumericVector pmf_min_range(int n, const std::vector<double>& pi,
                                  double c_lo, double c_hi) {
  if (!std::isfinite(c_lo) || !std::isfinite(c_hi))
    return Rcpp::NumericVector(0);

  const int lo = std::max(0, static_cast<int>(std::ceil(c_lo)));
  const int hi = std::min(n, static_cast<int>(std::floor(c_hi)));
  if (lo > hi) return Rcpp::NumericVector(0);

  const int len = hi - lo + 1;
  Rcpp::NumericVector out(len);
  // surv[i] = P(min >= lo + i),  surv[len] = P(min >= hi + 1)
  std::vector<double> surv(len + 1);
  for (int i = 0; i <= len; ++i)
    surv[i] = prob_min_geq(n, pi, static_cast<double>(lo + i));
  for (int i = 0; i < len; ++i)
    out[i] = std::max(0.0, surv[i] - surv[i + 1]);
  return out;
}

static double prob_range_lt_int(int n, const std::vector<double>& pi, int r) {
  if (r <= 0) return 0.0;
  if (r > n)  return 1.0;
  const int m = static_cast<int>(pi.size());
  if (m <= 1) return 1.0;   // single urn: range is always 0

  double result = 0.0;

  // First sum: h = 0..n-r+1,  a = h+r-1, b = h
  for (int h = 0; h <= n - r + 1; ++h)
    result += prob_joint(n, pi,
                         static_cast<double>(h + r - 1),
                         static_cast<double>(h));

  // Second sum (subtracted): h = 0..n-r,  a = h+r-1, b = h+1
  for (int h = 0; h <= n - r; ++h)
    result -= prob_joint(n, pi,
                         static_cast<double>(h + r - 1),
                         static_cast<double>(h + 1));

  return std::min(1.0, std::max(0.0, result));
}

// [[Rcpp::export]]
double prob_range_lt(int n, const std::vector<double>& pi, double r) {
  if (std::isnan(r))              return 0.0;
  if (r <= 0.0)                   return 0.0;
  if (r > static_cast<double>(n)) return 1.0;
  return prob_range_lt_int(n, pi, static_cast<int>(std::ceil(r)));
}

// [[Rcpp::export]]
double prob_range_leq(int n, const std::vector<double>& pi, double r) {
  if (std::isnan(r))               return 0.0;
  if (r < 0.0)                     return 0.0;
  if (r >= static_cast<double>(n)) return 1.0;
  return prob_range_lt_int(n, pi, static_cast<int>(std::floor(r)) + 1);
}

// [[Rcpp::export]]
double prob_range_geq(int n, const std::vector<double>& pi, double r) {
  return std::max(0.0, 1.0 - prob_range_lt(n, pi, r));
}

// [[Rcpp::export]]
double prob_range_gt(int n, const std::vector<double>& pi, double r) {
  return std::max(0.0, 1.0 - prob_range_leq(n, pi, r));
}

// Pr( range = r )   —   PMF; returns 0 for non-integer r
// P(range = r) = P(range < r+1) - P(range < r)
// [[Rcpp::export]]
double prob_range_eq(int n, const std::vector<double>& pi, double r) {
  if (!std::isfinite(r) || r != std::floor(r)) return 0.0;
  const int ri = static_cast<int>(r);
  if (ri < 0 || ri > n) return 0.0;
  double hi = prob_range_lt_int(n, pi, ri + 1);
  double lo = (ri > 0) ? prob_range_lt_int(n, pi, ri) : 0.0;
  return std::max(0.0, hi - lo);
}

// PMF over a range of integer r values [ceil(r_lo), floor(r_hi)]
// [[Rcpp::export]]
Rcpp::NumericVector pmf_range_range(int n, const std::vector<double>& pi,
                                    double r_lo, double r_hi) {
  if (!std::isfinite(r_lo) || !std::isfinite(r_hi))
    return Rcpp::NumericVector(0);

  const int lo = std::max(0, static_cast<int>(std::ceil(r_lo)));
  const int hi = std::min(n, static_cast<int>(std::floor(r_hi)));
  if (lo > hi) return Rcpp::NumericVector(0);

  const int len = hi - lo + 1;
  Rcpp::NumericVector out(len);

  // CDF anchor just below lo: P(range < lo)
  double prev = (lo > 0) ? prob_range_lt_int(n, pi, lo) : 0.0;
  for (int i = 0; i < len; ++i) {
    double curr = prob_range_lt_int(n, pi, lo + i + 1);
    out[i] = std::max(0.0, curr - prev);
    prev = curr;
  }
  return out;
}
