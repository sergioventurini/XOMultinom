// -*- mode: C++; c-indent-level: 2; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

//' CDF of the sum of the first \eqn{J} largest order statistics for a
//' multinomial distribution.
//'
//' This function calculates the cumulative distribution function (CDF) of the
//' sum of the first \eqn{J} largest order statistics under an equiprobable
//' multinomial distribution assumption.
//'
//' @param t A length-one numeric vector indicating the value to compute the CDF
//'   for.
//' @param n A length-one integer vector indicating the number of independent
//'   balls.
//' @param m A length-one integer vector indicating the number of independent
//'   urns/cells.
//' @param J A length-one integer vector indicating the number of largest order
//'   statistics.
//' @return A length-one numeric vector representing the probability of the
//'   sum of the \eqn{J} largest order statistic.
//' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
//' @seealso \code{\link{max_order_statistic}} for computing the
//'   CDF of the maximum.
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
//' highest_order_statistics_C(6, 10, 5, 2)
//'
// [[Rcpp::export]]
double highest_order_statistics_C(const double & td, int n, int m, int J) {
  int t = floor(td);
  double P = 0;

  if (J > m) {
    Rcpp::stop("J should be smaller or equal than m.");
  }
  if ((J == m) && (t < n)) {
    Rcpp::stop("total sum is every time equal to n.");
  }
  if ((t == 0) && (n != 0)) {
    P = 0;
    return P;
  }
  if ((t == 0) && (n == 0)) {
    P = 1;
    return P;
  }

  // first term if n_<1> <= t/J
  P = max_order_statistic_C(floor(t/J), n, m);

  if (J > 1) {
    for (int sum_depth = 1; sum_depth <= (J - 1); sum_depth++) {
      arma::vec rangeArg = arma::vec();
      int cur_depth = 1;
      // double TMP = recursive_sum_C(t, n, m, J, sum_depth, cur_depth, rangeArg);
      // printf("recursive_sum_C(%.1d, %.0d, %.0d, %.0d, %.0d, %.0d, rangeArg) = %.7f\n", t, n, m, J, sum_depth, cur_depth, TMP);
      P = P + recursive_sum_C(t, n, m, J, sum_depth, cur_depth, rangeArg);
      R_CheckUserInterrupt();
    }
  }

  return P;
}
