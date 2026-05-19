#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

/*
//' CDF of the sum of the first \eqn{J} largest order statistics
//'
//' Computes the cumulative distribution function of the sum of the first
//' \eqn{J} largest order statistics for an equiprobable multinomial
//' distribution.
//'
//' @param td Numeric value at which to evaluate the CDF. Internally, this is
//'   converted to \code{floor(td)}.
//' @param n Integer number of balls/trials.
//' @param m Integer number of urns/cells.
//' @param J Integer number of largest order statistics to sum.
//'
//' @return Numeric value giving the CDF probability.
//'
//' @references
//' Bonetti, M., Cirillo, P., Ogay, A. (2019). Computing the exact
//' distributions of some functions of the ordered multinomial counts:
//' maximum, minimum, range and sums of order statistics. Royal Society
//' Open Science, 6, 190198. \doi{10.1098/rsos.190198}
//'
//' @examples
//' highest_order_statistics(6, 10, 5, 2)
//'
//' @keywords internal
//'
*/
// [[Rcpp::export]]
double highest_order_statistics(const double & td, int n, int m, int J) {
  int t = (int)floor(td);
  double P = 0.0;

  if (J > m) {
    Rcpp::stop("J must be less than or equal to m.");
  }
  if ((J == m) && (t < n)) {
    Rcpp::stop("total sum always equals n.");
  }
  if ((t == 0) && (n != 0)) {
    return 0.0;
  }
  if ((t == 0) && (n == 0)) {
    return 1.0;
  }

  P = max_order_statistic((double)floor((double)t / J), n, m);

  if (J > 1) {
    std::vector<int> rangeArg;
    rangeArg.reserve(J - 1);   // maximum number of elements pushed at any time

    for (int sum_depth = 1; sum_depth <= (J - 1); sum_depth++) {
      rangeArg.clear();
      P += recursive_sum_impl(t, n, m, J, sum_depth,
                              /* cur_depth = */ 1,
                              rangeArg,
                              /* partial_sum = */ 0);
      R_CheckUserInterrupt();
    }
  }

  return P;
}
