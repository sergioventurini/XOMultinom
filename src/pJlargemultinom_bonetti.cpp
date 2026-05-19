#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

/*
//' CDF of the sum of \eqn{J} largest order statistics for a multinomial
//' distribution (equiprobable case)
//'
//' Computes the cumulative distribution function of the sum of \eqn{J}
//' largest order statistics for a multinomial distribution with
//' equal cell probabilities.
//'
//' @param x Numeric vector of values at which to evaluate the CDF.
//' @param size Integer number of trials.
//' @param prob Numeric vector of category probabilities. Must sum to 1.
//' @param J Integer number of order statistics.
//' @param logd Logical; if \code{TRUE}, returns log-probabilities.
//' @param verbose Logical; if \code{TRUE}, displays a progress bar.
//'
//' @return Numeric vector of the same length as \code{x} containing
//'   \eqn{P(S_J \le x)} (or log-probabilities if \code{logd = TRUE}).
//'
//' @examples
//' pJlargemultinom_bonetti(0:3, size = 10, prob = rep(1/5, 5),
//'                         J = 3, logd = FALSE, verbose = FALSE)
//'
//' @keywords internal
//'
*/
// [[Rcpp::export]]
Rcpp::NumericVector pJlargemultinom_bonetti(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob, const int& J,
  const bool& logd, const bool& verbose) {
  int xlen = x.size();
  Progress prog(xlen, verbose);  // verbose = false silences the bar

  int m = prob.size();

  Rcpp::NumericVector r(xlen);
  for (int k = 0; k < xlen; k++) {
    if (Progress::check_abort()) return Rcpp::NumericVector(0);
    prog.increment();
    r(k) = highest_order_statistics(x(k), size, m, J);

    R_CheckUserInterrupt();
  }

  if (!logd) {
    return r;
  } else {
    return Rcpp::log(r);
  }
}
