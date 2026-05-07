#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

/*
//' Internal utility for smallest order statistic computation
//'
//' Recursive helper used to compute probabilities related to the
//' smallest order statistic in a multinomial distribution.
//'
//' @param t_max_d Numeric upper bound (internally coerced to integer).
//' @param n Integer number of balls/trials.
//' @param m Integer number of urns/cells.
//' @param t Integer threshold for the minimum.
//'
//' @return Numeric value.
//'
//' @keywords internal
//'
*/
// [[Rcpp::export]]
double max_for_min(const double & t_max_d, int n, int m, int t) {
  int t_max = (int)t_max_d;

  if (t_max < t) {
    return (n == 0 && m == 0) ? 1.0 : 0.0;
  }
  if (n == 0 && m == 0) return 1.0;
  if (t_max == 1) {
    // All cells have exactly 0 or 1 ball; minimum >= t requires m == n
    return (m == n)
             ? std::exp(lgamma(m + 1) - lgamma(m - n + 1) - (double)n * std::log((double)m))
             : 0.0;
  }
  if (n == 0 && m != 0) return 0.0;

  static std::unordered_map<std::tuple<int,int,int,int>, double, TupleHash4> cache;

  auto key = std::make_tuple(t_max, n, m, t);
  {
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;
  }

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
      aux += coef * max_for_min((double)(t_max - 1), n - t_max * q, m - q, t);
      R_CheckUserInterrupt();
    }
  }

  cache[key] = aux;
  return aux;
}
