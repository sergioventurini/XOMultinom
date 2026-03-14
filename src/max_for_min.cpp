#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

//' Utility function for computing the distribution of the smallest order
//' statistic.
//'
//' This is an auxiliary function to the distribution of the smallest order
//' statistic.
//'
//' @param t_max A length-one numeric vector indicating the value to compute the
//'   survival function for.
//' @param n A length-one integer vector indicating the number of independent
//'   balls.
//' @param m A length-one integer vector indicating the number of independent
//'   urns/cells.
//' @param t A length-one numeric vector indicating the value to compute for.
//' @return A length-one numeric vector.
//' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
//' @seealso \code{\link{highest_order_statistics}} for computing the
//'   CDF of the sum of the first \eqn{J} largest order statistics.
//' @seealso \code{\link{smallest_order_value}} for computing the CDF
//'   of the smallest order statistic.
//' @seealso \code{\link{range_probability}} for computing the CDF
//'   of the range.
//' @references
//'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
//'   distributions of some functions of the ordered multinomial counts:
//'   maximum, minimum, range and sums of order statistics", Royal Society
//'   Open Science, 6: 190198, <http://dx.doi.org/10.1098/rsos.190198>.
//' @examples
//' smallest_order_value(1, 10, 5) # P(N_(1) <= 1; n = 10, m = 5)
//'
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
