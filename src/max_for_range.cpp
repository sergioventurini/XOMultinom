#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

double max_for_range_impl(double t_max_d, int n, int m,
                          int orig_t_max, int orig_n, int orig_m, int t) {
  int t_max = (int)t_max_d;

  if (t_max == orig_t_max && n == orig_n && m == orig_m) {
    return 0.0;
  }

  if (orig_t_max + 1 - t_max > t) {
    return (n == 0 && m == 0) ? 1.0 : 0.0;
  }

  if (n == 0 && m == 0) return 1.0;
  if (t_max == 1) {
    return (m == n)
             ? std::exp(lgamma(m + 1) - lgamma(m - n + 1) - (double)n * std::log((double)m))
             : 0.0;
  }
  if (n == 0 && m != 0) return 0.0;

  double aux = 0.0;
  const double common_term = lgamma(n + 1) + lgamma(m + 1) - (double)n * std::log((double)m);
  const int LowSum = (int)std::max(0, n - t_max * m + m);
  const int UpSum  = (int)floor((double)n / (double)t_max);

  if (LowSum <= UpSum) {
    for (int q = LowSum; q <= UpSum; q++) {
      const double summ_term = -q * lgamma(t_max + 1)
                               - lgamma(q + 1)
                               - lgamma(m - q + 1)
                               - lgamma(n - t_max * q + 1);
      const double summ_nom  = (m != q)
                                 ? (double)(n - t_max * q) * std::log((double)(m - q))
                                 : 0.0;
      const double coef = std::exp(common_term + summ_term + summ_nom);
      // Pass orig_* unchanged: `prev` is fixed throughout the entire recursion tree
      aux += coef * max_for_range_impl((double)(t_max - 1), n - t_max * q, m - q,
                                       orig_t_max, orig_n, orig_m, t);
      R_CheckUserInterrupt();
    }
  }

  return aux;
}

/*
//' Internal utility for range distribution computation
//'
//' Wrapper around the internal recursive implementation used to compute
//' probabilities related to the range of multinomial cell counts.
//'
//' @param t_max Numeric upper bound.
//' @param n Integer number of balls/trials.
//' @param m Integer number of urns/cells.
//' @param prev Numeric vector of length 3 containing previous recursion values.
//' @param t Integer threshold value.
//'
//' @return Numeric probability value.
//'
//' @keywords internal
//'
*/
// [[Rcpp::export]]
double max_for_range(const double & t_max, int n, int m,
  Rcpp::NumericVector prev, int t) {
  return max_for_range_impl(t_max, n, m,
                            (int)prev[0], (int)prev[1], (int)prev[2],
                            t);
}
