// -*- mode: C++; c-indent-level: 2; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <R_ext/Utils.h>

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
//' smallest_order_value_C(1, 10, 5) # P(N_(1) <= 1; n = 10, m = 5)
//'
// [[Rcpp::export]]
double max_for_min_C(const double & t_max, int n, int m, int t) {
  double aux = 0;
  if (t_max < t) {
    if (n == 0 && m == 0) {
      aux = 1;
      return aux;
    } else {
      aux = 0;
      return aux;
    }
  } else {
    if (n == 0 && m == 0) {
      aux = 1;
      return aux;
    }
    if (t_max == 1) {
      if (m == n) {
        aux = exp(lgamma(m + 1) - lgamma(m - n + 1) - n*log(m)); // explicit calculation of P(n<1> <= 1);
        return aux;
      } else {
        aux = 0;
        return aux;
      }
    }
    if (n == 0 && m != 0) {
      aux = 0;
      return aux;
    }
    double common_term = lgamma(n + 1) + lgamma(m + 1) - n*log(m);
    int LowSum = fmax(0, n - t_max*m + m);
    int UpSum = floor(n/t_max);
    if (LowSum <= UpSum) {
      for (int q = LowSum; q <= UpSum; q++) {
        double summ_term = -q*lgamma(t_max + 1) - lgamma(q + 1) - lgamma(m - q + 1) - lgamma(n - t_max*q + 1);
        double summ_term_nominator = 0;
        if (m != q) {
          summ_term_nominator = (n - t_max*q)*log(m - q);
        }
        double coef = exp(common_term + summ_term + summ_term_nominator);
        double temp = max_for_min_C(t_max - 1, n - t_max*q, m - q, t);
        aux = aux + coef*temp;
        R_CheckUserInterrupt();
      }
    }
  }

  return aux;
}
