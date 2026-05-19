#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

/*
//' PMF of the sum of \eqn{J} largest order statistics for a multinomial
//' distribution (equiprobable case)
//'
//' Computes the probability mass function of the sum of the \eqn{J} largest
//' order statistics for a multinomial random vector with total count `size`
//' and equal category probabilities `prob`, evaluated at each value in `x`.
//'
//' @param x Numeric vector of values at which to evaluate the probability mass
//'   function.
//' @param size Integer total count of the multinomial distribution.
//' @param prob Numeric vector of category probabilities.
//' @param J Integer number of order statistics.
//' @param logd Logical; if `TRUE`, probabilities are returned on the log scale.
//' @param verbose Logical; if `TRUE`, displays a progress bar.
//'
//' @return A numeric vector containing the probabilities, or log-probabilities
//'   when `logd = TRUE`.
//'
//' @keywords internal
//'
*/
// [[Rcpp::export]]
Rcpp::NumericVector dJlargemultinom_bonetti(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob, const int& J,
  const bool& logd, const bool& verbose) {
  int xlen = x.size();
  Progress prog(xlen, verbose);  // verbose = false silences the bar

  int m = prob.size();

  Rcpp::NumericVector r(xlen);
  for (int k = 0; k < xlen; k++) {
    if (Progress::check_abort()) return Rcpp::NumericVector(0);
    prog.increment();
    r(k) = equi_prob_highest_eq(size, m, x(k), J);

    R_CheckUserInterrupt();
  }

  if (!logd) {
    return r;
  } else {
    return Rcpp::log(r);
  }
}
