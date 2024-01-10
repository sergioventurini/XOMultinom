// -*- mode: C++; c-indent-level: 2; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

//' CDF of the range for an equiprobable multinomial distribution.
//'
//' This function calculates the cumulative distribution function (CDF)
//' of the range under the assumption of an equiprobable multinomial
//' distribution.
//'
//' @param td A length-one numeric vector indicating the value to compute the
//'   survival function for.
//' @param n A length-one integer vector indicating the number of independent
//'   balls.
//' @param m A length-one integer vector indicating the number of independent
//'   urns/cells.
//' @return A length-one numeric vector representing the probability of the
//'   smallest order statistic.
//' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
//' @seealso \code{\link{highest_order_statistics}} for computing the
//'   CDF of the sum of the first \eqn{J} largest order statistics.
//' @seealso \code{\link{smallest_order_value}} for computing the CDF
//'   of the smallest order statistic.
//' @seealso \code{\link{highest_order_statistics}} for computing the sum of
//'   the first \eqn{J} largest order statistics.
//' @references
//'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
//'   distributions of some functions of the ordered multinomial counts:
//'   maximum, minimum, range and sums of order statistics", Royal Society
//'   Open Science, 6: 190198, <http://dx.doi.org/10.1098/rsos.190198>.
//' @examples
//' range_probability_C(1, 7, 10)
//'
// [[Rcpp::export]]
double range_probability_C(const double & td, int n, int m) {
  int t = floor(td);
  double P = 0;
  double aux = 0;

  if (t > n) {
    P = 1;
    return P;
  }
  P = max_order_statistic_C(t, n, m);
  arma::vec prev = arma::vec(3);
  prev(0) = t;
  prev(1) = n;
  prev(2) = m;
  for (int t_max = (t + 1); t <= n; t++) {
    aux = max_for_range_C(t_max, n, m, prev, t);
    P = P + aux;
    prev(0) = t_max;
    prev(1) = n;
    prev(2) = m;
    R_CheckUserInterrupt();
  }

  return P;
}
