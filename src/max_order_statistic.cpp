// -*- mode: C++; c-indent-level: 2; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <R_ext/Utils.h>

// Note: RcppExport is an alias for extern "C"

//' CDF of the highest order statistic for an equiprobable multinomial
//' distribution.
//'
//' This function calculates the cumulative distribution function (CDF) of the
//' highest order statistic (i.e. the maximum) under an equiprobable multinomial
//' distribution assumption.
//'
//' @param td A length-one numeric vector indicating the value to compute the CDF
//'   for.
//' @param n A length-one integer vector indicating the number of independent
//'   balls.
//' @param m A length-one integer vector indicating the number of independent
//'   urns/cells.
//' @return A length-one numeric vector representing the probability of the
//'   highest order statistic.
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
//' max_order_statistic_C(3, 10, 5) # P(N_(m) <= 3; n = 10, m = 5)
//'
// [[Rcpp::export]]
double max_order_statistic_C(const double & td, int n, int m) {
  int t = floor(td);
  double P = 0;

  if ((t == 0) && (n != 0)) {
    P = 0;
    return P;
  }
  if ((t == 0) && (n == 0)) {
    P = 1;
    return P;
  }
  if (t >= n) {
    P = 1;
    return P;
  }
  if ((n == 0) && (m != 0)) {
    P = 1;
    return P;
  }
  if ((n == 0) && (m == 0)) {
    P = 1;
    return P;
  }

  double common_term = lgamma(n + 1) + lgamma(m + 1) - n*log(m);
  if (t == 1) {
    if (m >= n) {
      P = exp(lgamma(m + 1) - lgamma(m - n + 1) - n*log(m));
    } else {
      P = 0;
    }
  } else {
    int LowSum = fmax(0, n - t*m + m);
    int UpSum = floor(n/t);
    if (LowSum <= UpSum) {
      for (int q = LowSum; q <= UpSum; q++) {
        double summ_term = -q*lgamma(t + 1) - lgamma(q + 1) - lgamma(m - q + 1) - lgamma(n - t*q + 1);
        double summ_term_nominator = 0;
        if (m != q) {
          summ_term_nominator = (n - t*q)*log(m - q);
        }
        double coef = exp(common_term + summ_term + summ_term_nominator);
        double temp = max_order_statistic_C(t - 1, n - t*q, m - q);
        P = P + coef*temp;
        R_CheckUserInterrupt();
      }
    }
  }

  return P;
}
