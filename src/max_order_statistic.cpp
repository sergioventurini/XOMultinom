#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

/*
//' CDF of the maximum multinomial cell count
//'
//' Computes the cumulative distribution function of the largest order
//' statistic, that is, the maximum cell count, for an equiprobable
//' multinomial distribution.
//'
//' @param td Numeric value at which to evaluate the CDF. Internally, this is
//'   converted to \code{floor(td)}.
//' @param n Integer number of balls/trials.
//' @param m Integer number of urns/cells.
//'
//' @return Numeric value giving \eqn{P(N_{(m)} \le t)}, where
//'   \eqn{t = \lfloor td \rfloor}.
//'
//' @references
//' Bonetti, M., Cirillo, P., Ogay, A. (2019). Computing the exact
//' distributions of some functions of the ordered multinomial counts:
//' maximum, minimum, range and sums of order statistics. Royal Society
//' Open Science, 6, 190198. \doi{10.1098/rsos.190198}
//'
//' @examples
//' max_order_statistic(3, 10, 5) # P(N_(m) <= 3; n = 10, m = 5)
//'
//' @keywords internal
//'
*/
// [[Rcpp::export]]
double max_order_statistic(const double & td, int n, int m) {
  int t = (int)floor(td);

  if (t == 0 && n != 0) return 0.0;
  if (t >= n) return 1.0;

  // --- Memoization ---
  // Static cache persists across R calls, so subproblems are reused when the
  // user evaluates the CDF at multiple points or when called from
  // highest_order_statistics / range_probability.
  static std::unordered_map<std::tuple<int,int,int>, double, TupleHash3> cache;

  auto key = std::make_tuple(t, n, m);
  {
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;
  }

  // --- Recursive computation ---
  double P = 0.0;
  const double common_term = lgamma(n + 1) + lgamma(m + 1) - (double)n * std::log((double)m);

  if (t == 1) {
    // P(all cells <= 1 | n, m): requires m >= n (pigeon-hole)
    P = (m >= n)
          ? std::exp(lgamma(m + 1) - lgamma(m - n + 1) - (double)n * std::log((double)m))
          : 0.0;
  } else {
    const int LowSum = (int)std::max(0, n - t * m + m);
    const int UpSum  = (int)floor((double)n / (double)t);

    for (int q = LowSum; q <= UpSum; q++) {
      const double summ_term = -q * lgamma(t + 1)
                               - lgamma(q + 1)
                               - lgamma(m - q + 1)
                               - lgamma(n - t * q + 1);
      // Guard log(0): when m == q all remaining cells have exactly t balls
      const double summ_nom  = (m != q)
                                 ? (double)(n - t * q) * std::log((double)(m - q))
                                 : 0.0;
      P += std::exp(common_term + summ_term + summ_nom)
           * max_order_statistic((double)(t - 1), n - t * q, m - q);
      R_CheckUserInterrupt();
    }
  }

  cache[key] = P;
  return P;
}
