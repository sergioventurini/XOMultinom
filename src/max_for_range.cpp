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
      aux += coef * max_for_range_impl((double)(t_max - 1), m - t_max * q, m - q,
                                       orig_t_max, orig_n, orig_m, t);
      R_CheckUserInterrupt();
    }
  }

  return aux;
}

//' Utility function for computing the distribution of the range.
//'
//' This is an auxiliary function to the distribution of the range.
//'
//' @param t_max A length-one numeric vector indicating the value to compute the
//'   survival function for.
//' @param n A length-one integer vector indicating the number of independent
//'   balls.
//' @param m A length-one integer vector indicating the number of independent
//'   urns/cells.
//' @param prev Integer vector containing the previous iteration values.
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
//' range_probability(1, 7, 10)
//'
// [[Rcpp::export]]
double max_for_range(const double & t_max, int n, int m, arma::vec prev, int t) {
  return max_for_range_impl(t_max, n, m,
                            (int)prev(0), (int)prev(1), (int)prev(2),
                            t);
}
