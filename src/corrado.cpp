#include "XOMultinom.h"

// [[Rcpp::depends("RcppArmadillo")]]

// Note: RcppExport is an alias for extern "C"

static inline double log_binom_pmf(int N, double p, int j) {
  if (j < 0 || j > N) return -std::numeric_limits<double>::infinity();
  if (p <= 0.0)        return (j == 0) ? 0.0 : -std::numeric_limits<double>::infinity();
  if (p >= 1.0)        return (j == N) ? 0.0 : -std::numeric_limits<double>::infinity();
  return std::lgamma(N + 1.0)
       - std::lgamma(j + 1.0)
       - std::lgamma(N - j + 1.0)
       + j  * std::log(p)
       + (N - j) * std::log(1.0 - p);
}

// fill_binom_range --> fills v[lo..hi] with Binom(N, p, j) for j in [lo, hi]
static void fill_binom_range(std::vector<double>& v,
                             int N, double p, int lo, int hi) {
  if (lo > hi) return;

  const double q  = 1.0 - p;
  const double pq = (q > 1e-300) ? p / q : 0.0;
  const double qp = (p > 1e-300) ? q / p : 0.0;

  if (p > 1.0 - 1e-14) {
    for (int j = lo; j <= hi; ++j) v[j] = (j == N) ? 1.0 : 0.0;
    return;
  }
  if (p < 1e-300) {
    for (int j = lo; j <= hi; ++j) v[j] = (j == 0) ? 1.0 : 0.0;
    return;
  }

  if (lo == 0) {
    double f = std::pow(q, static_cast<double>(N));
    if (f > 0.0) {
      v[0] = f;
      for (int j = 1; j <= hi; ++j) {
        f    *= pq * (N - j + 1.0) / j;
        v[j]  = f;
      }
      return;
    }
  }

  // ── General path (lo > 0, or lo = 0 with pow underflow) ─────────────────
  // Seed at the binomial mode clipped to [lo, hi]; always a normal double.

  const int mode  = static_cast<int>(std::floor(static_cast<double>(N) * p));
  const int jstar = std::max(lo, std::min(hi, mode));

  double f  = std::exp(log_binom_pmf(N, p, jstar));
  v[jstar]  = f;

  // Walk UP: Binom(N,p,j+1) = Binom(N,p,j) * (p/q) * (N-j)/(j+1)
  double fu = f;
  for (int j = jstar; j < hi; ++j) {
    fu       *= pq * (N - j) / (j + 1.0);
    v[j + 1]  = fu;
  }

  // Walk DOWN: Binom(N,p,j-1) = Binom(N,p,j) * (q/p) * j/(N-j+1)
  double fd = f;
  for (int j = jstar; j > lo; --j) {
    fd       *= qp * j / (N - j + 1.0);
    v[j - 1]  = fd;
  }
}

// vstep --> convolve weight vector v[] with Binom(n-i, p, d) for d in [b, a]
static void vstep(std::vector<double>& v, std::vector<double>& w,
                  int n, double p, int b, int a) {
  const double q  = 1.0 - p;
  const double pq = (q > 1e-300) ? p / q : 0.0;

  std::fill(w.begin(), w.end(), 0.0);

  // Edge case: p = 1  =>  Binom(n-i, 1, d) puts all mass at d = n-i
  if (p > 1.0 - 1e-14) {
    for (int i = 0; i <= n; ++i) {
      if (v[i] == 0.0) continue;
      const int d = n - i;
      if (d >= b && d <= a) w[n] += v[i];
    }
    std::swap(v, w);
    return;
  }

  // precompute g[i] = Binom(n-i, p, b) for all i

  std::vector<double> g(n + 1, 0.0);

  if (b == 0) {
    // fast path (b = 0, max statistic)
    g[n] = 1.0;
    for (int i = n; i > 0; --i) {
      g[i - 1] = g[i] * q;
      if (g[i - 1] == 0.0) break;   // underflow: tail is zero, stop early
    }

  } else {
    // ── General path (b > 0, min and range statistics) ───────────────────────
    // Peak of g[i] as a function of i: i* = n - b/p  (clipped to [0, n])
    const int istar = std::max(0, std::min(n,
          n - static_cast<int>(std::ceil(static_cast<double>(b) / p))));

    // seed at istar in log-space (always a normal double, immune to underflow)
    const int max_d_star = n - istar;
    if (max_d_star >= b) {
      g[istar] = std::exp(log_binom_pmf(max_d_star, p, b));

      // walk UP:   g[i+1] = g[i] * (n-i-b) / ((n-i)*q)
      for (int i = istar; i < n; ++i) {
        if (g[i] == 0.0) break;
        const int mdi = n - i;
        if (mdi <= b) break;
        g[i + 1] = g[i] * static_cast<double>(mdi - b)
                         / (static_cast<double>(mdi) * q);
      }

      // Walk DOWN: g[i-1] = g[i] * (n-i+1)*q / (n-i+1-b)
      for (int i = istar; i > 0; --i) {
        if (g[i] == 0.0) break;
        const int mdi1 = n - i + 1;
        if (mdi1 <= b) break;
        g[i - 1] = g[i] * static_cast<double>(mdi1) * q
                         / static_cast<double>(mdi1 - b);
      }
    }
  }

  // main convolution loop
  for (int i = 0; i <= n; ++i) {
    if ((i & 1023) == 0) R_CheckUserInterrupt();   // throttled: every 1024 rows
    if (v[i] == 0.0 || g[i] == 0.0) continue;
    const int max_d = n - i;
    // g[i] > 0 implies max_d >= b (ensured by the precomputation above)

    double f         = g[i];
    const int  d_hi  = std::min(a, max_d);
    const double vi  = v[i];
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
  if (!max_convert(c, n, ci)) return 0.0;
  if (ci >= n) return 1.0;
  if (ci <  0) return 0.0;

  const int m = static_cast<int>(pi.size());

  if (static_cast<long long>(ci) * m < static_cast<long long>(n)) return 0.0;

  std::vector<double> suf(m + 1, 0.0);
  for (int k = m - 1; k >= 0; --k) suf[k] = suf[k + 1] + pi[k];

  std::vector<double> v(n + 1, 0.0), w(n + 1);

  // Initialise v[j] = Binom(n, pi[0], j) for j in [0, ci].
  // fill_binom_range seeds at the mode so pow(q,n) is never computed.
  {
    const double p = pi[0];
    if (p > 1.0 - 1e-14) { if (n <= ci) v[n] = 1.0; }
    else                  { fill_binom_range(v, n, p, 0, std::min(ci, n)); }
  }

  for (int k = 1; k <= m - 2; ++k)
    vstep(v, w, n, pi[k] / suf[k], 0, ci);

  double result = 0.0;
  for (int j = std::max(0, n - ci); j <= n; ++j) result += v[j];
  return std::min(1.0, result);
}

// [[Rcpp::export]]
double prob_min_geq(int n, const std::vector<double>& pi, double c) {
  int ci;
  if (!min_convert(c, n, ci)) return 1.0;
  if (ci <= 0)                return 1.0;
  if (static_cast<long long>(ci) * pi.size() > static_cast<long long>(n)) return 0.0;

  const int m = static_cast<int>(pi.size());

  std::vector<double> suf(m + 1, 0.0);
  for (int k = m - 1; k >= 0; --k) suf[k] = suf[k + 1] + pi[k];

  std::vector<double> v(n + 1, 0.0), w(n + 1);

  // Initialise v[j] = Binom(n, pi[0], j) for j in [ci, n].
  {
    const double p = pi[0];
    if (p > 1.0 - 1e-14) { if (n >= ci) v[n] = 1.0; }
    else                  { fill_binom_range(v, n, p, ci, n); }
  }

  for (int k = 1; k <= m - 2; ++k)
    vstep(v, w, n, pi[k] / suf[k], ci, n);

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
  if (c < 0.0)       return 0.0;
  if (c >= static_cast<double>(n)) return 1.0;
  return std::max(0.0, 1.0 - prob_min_geq(n, pi, std::floor(c) + 1.0));
}

// [[Rcpp::export]]
double prob_joint(int n, const std::vector<double>& pi, double a, double b) {
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

  // Initialise v[j] = Binom(n, pi[0], j) for j in [lo, hi].
  {
    const double p  = pi[0];
    const int lo = std::max(0, bi), hi = std::min(ai, n);
    if (p > 1.0 - 1e-14) { if (n >= lo && n <= hi) v[n] = 1.0; }
    else                  { fill_binom_range(v, n, p, lo, hi); }
  }

  for (int k = 1; k <= m - 2; ++k)
    vstep(v, w, n, pi[k] / suf[k], bi, ai);

  double result = 0.0;
  for (int j = std::max(0, n - ai); j <= std::min(n, n - bi); ++j)
    result += v[j];
  return std::min(1.0, result);
}

// PMF: P(max n_k = c)
// [[Rcpp::export]]
double prob_max_eq(int n, const std::vector<double>& pi, double c) {
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
// [[Rcpp::export]]
Rcpp::NumericVector pmf_max_range(int n, const std::vector<double>& pi,
                                  double c_lo, double c_hi) {
  if (!std::isfinite(c_lo) || !std::isfinite(c_hi))
    return Rcpp::NumericVector(0);

  const int lo = std::max(0, static_cast<int>(std::ceil(c_lo)));
  const int hi = std::min(n, static_cast<int>(std::floor(c_hi)));
  if (lo > hi) return Rcpp::NumericVector(0);

  const int len = hi - lo + 1;
  Rcpp::NumericVector out(len);
  double prev = (lo > 0) ? prob_max_leq(n, pi, static_cast<double>(lo - 1)) : 0.0;
  for (int i = 0; i < len; ++i) {
    R_CheckUserInterrupt();
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
  std::vector<double> surv(len + 1);
  for (int i = 0; i <= len; ++i) {
    R_CheckUserInterrupt();
    surv[i] = prob_min_geq(n, pi, static_cast<double>(lo + i));
  }
  for (int i = 0; i < len; ++i)
    out[i] = std::max(0.0, surv[i] - surv[i + 1]);
  return out;
}

static double prob_range_lt_int(int n, const std::vector<double>& pi, int r) {
  if (r <= 0) return 0.0;
  if (r > n)  return 1.0;
  const int m = static_cast<int>(pi.size());
  if (m <= 1) return 1.0;

  double result = 0.0;

  for (int h = 0; h <= n - r + 1; ++h) {
    R_CheckUserInterrupt();
    result += prob_joint(n, pi,
                         static_cast<double>(h + r - 1),
                         static_cast<double>(h));
  }

  for (int h = 0; h <= n - r; ++h) {
    R_CheckUserInterrupt();
    result -= prob_joint(n, pi,
                         static_cast<double>(h + r - 1),
                         static_cast<double>(h + 1));
  }

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

// [[Rcpp::export]]
double prob_range_eq(int n, const std::vector<double>& pi, double r) {
  if (!std::isfinite(r) || r != std::floor(r)) return 0.0;
  const int ri = static_cast<int>(r);
  if (ri < 0 || ri > n) return 0.0;
  double hi = prob_range_lt_int(n, pi, ri + 1);
  double lo = (ri > 0) ? prob_range_lt_int(n, pi, ri) : 0.0;
  return std::max(0.0, hi - lo);
}

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

  double prev = (lo > 0) ? prob_range_lt_int(n, pi, lo) : 0.0;
  for (int i = 0; i < len; ++i) {
    R_CheckUserInterrupt();
    double curr = prob_range_lt_int(n, pi, lo + i + 1);
    out[i] = std::max(0.0, curr - prev);
    prev = curr;
  }
  return out;
}
