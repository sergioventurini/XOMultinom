// -*- mode: C++; c-indent-level: 2; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

//' Utility function for computing recursively a sum.
//'
//' This is an auxiliary function to compute a sum recursively.
//'
//' @param td A length-one numeric vector indicating the value to compute the
//'   function for.
//' @param n A length-one integer vector indicating the number of independent
//'   balls.
//' @param m A length-one integer vector indicating the number of independent
//'   urns/cells.
//' @param J A length-one integer vector indicating the number of largest order
//'   statistics.
//' @param sum_depth A length-one integer vector indicating the summation depth.
//' @param cur_depth A length-one integer vector indicating the current
//'   summation depth.
//' @param rangeArg Integer vector indicating the range of the summation.
//' @return A length-one numeric vector representing the summation value.
//' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
//' @seealso \code{\link{highest_order_statistics}} for computing the CDF
//'   of the sum for the first \eqn{J} largest order statistics for an
//'   equiprobable multinomial distribution.
//' @references
//'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
//'   distributions of some functions of the ordered multinomial counts:
//'   maximum, minimum, range and sums of order statistics", Royal Society
//'   Open Science, 6: 190198, <http://dx.doi.org/10.1098/rsos.190198>.
//' @examples
//' highest_order_statistics(6, 10, 5, 2)
//'
// [[Rcpp::export]]
double recursive_sum_C(const double & td, int n, int m, int J, int sum_depth, int cur_depth, arma::vec rangeArg) {
  int t = floor(td);

  double S = 0;
  arma::vec cur_range = arma::vec(2);

  if (cur_depth <= sum_depth) { // either increment summation depth or calculate the term
    if (cur_depth == 1) {
      cur_range(0) = floor(t/J + 1);
      cur_range(1) = t - sum_depth + 1;
    } else {
      cur_range(0) = floor((t - arma::sum(rangeArg(arma::span(0, cur_depth - 2))))/(J - cur_depth + 1) + 1);
      cur_range(1) = fmin(rangeArg(cur_depth - 2), t - arma::sum(rangeArg(arma::span(0, cur_depth - 2))));
    }
    for (int r = cur_range(0); r <= cur_range(1); r++) {
      rangeArg.resize(cur_depth);
      rangeArg(cur_depth - 1) = r;
      S = S + recursive_sum_C(t, n, m, J, sum_depth, cur_depth + 1, rangeArg);
      R_CheckUserInterrupt();
    }
  } else {
    double prob_arg = floor((t - arma::sum(rangeArg))/(J - sum_depth));
    double temp_p = max_order_statistic_C(prob_arg, n - arma::sum(rangeArg), m - sum_depth);
    double common_term = lgamma(n + 1) + lgamma(m + 1) - n*log(m);
    double coef = (n - arma::sum(rangeArg))*log(m - sum_depth) - lgamma(m - sum_depth + 1) - lgamma(n - arma::sum(rangeArg) + 1);
    for (int k = 0; k < rangeArg.n_elem; k++) {
      coef = coef - lgamma(rangeArg(k) + 1);
      R_CheckUserInterrupt();
    }
    arma::vec equal_statistics = arma::unique(rangeArg);
    for (int k = 0; k < equal_statistics.n_elem; k++) {
      arma::uvec temp_idx = arma::find(rangeArg == equal_statistics(k));
      int temp = temp_idx.n_elem;
      coef = coef - lgamma(temp + 1);
      R_CheckUserInterrupt();
    }
    S = temp_p*exp(common_term + coef);
  }

  return S;
}
