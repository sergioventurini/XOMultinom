#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

/*
//' CDF of the range of multinomial cell counts
//'
//' Computes the cumulative distribution function of the range, that is,
//' the difference between the maximum and minimum cell counts, for an
//' equiprobable multinomial distribution.
//'
//' @param td Numeric value at which to evaluate the CDF. Internally, this is
//'   converted to \code{floor(td)}.
//' @param n Integer number of balls/trials.
//' @param m Integer number of urns/cells.
//'
//' @return Numeric value giving \eqn{P(N_{(m)} - N_{(1)} \le t)}, where
//'   \eqn{t = \lfloor td \rfloor}.
//'
//' @references
//' Bonetti, M., Cirillo, P., Ogay, A. (2019). Computing the exact
//' distributions of some functions of the ordered multinomial counts:
//' maximum, minimum, range and sums of order statistics. Royal Society
//' Open Science, 6, 190198. \doi{10.1098/rsos.190198}
//'
//' @examples
//' range_probability(1, 7, 10)
//'
//' @keywords internal
//'
*/
// [[Rcpp::export]]
double range_probability(const double & td, int n, int m) {
  int t = (int)floor(td);
  double P = 0.0;

  if (t > n) {
    return 1.0;
  }

  P = max_order_statistic((double)t, n, m);

  int prev_t_max = t;   // initial "previous" t_max before the loop starts

  for (int t_max = t + 1; t_max <= n; t_max++) {
    P += max_for_range_impl((double)t_max, n, m,
                            prev_t_max, n, m,
                            t);
    prev_t_max = t_max;
    R_CheckUserInterrupt();
  }

  return P;
}
