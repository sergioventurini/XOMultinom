#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

/*
//' CDF of the minimum for a multinomial distribution
//'
//' Computes the cumulative distribution function of the minimum
//' \eqn{M = \min(X_1, \ldots, X_k)} for a multinomial distribution
//' with equal cell probabilities.
//'
//' @param x Numeric vector of values at which to evaluate the CDF.
//' @param size Integer number of trials.
//' @param prob Numeric vector of category probabilities. Must sum to 1.
//' @param logd Logical; if \code{TRUE}, returns log-probabilities.
//' @param verbose Logical; if \code{TRUE}, displays a progress bar.
//'
//' @return Numeric vector of the same length as \code{x} containing
//'   \eqn{P(M \le x)} (or log-probabilities if \code{logd = TRUE}).
//'
//' @examples
//' pminmultinom_bonetti(0:3, size = 10, prob = rep(1/5, 5),
//'                      logd = FALSE, verbose = FALSE)
//'
//' @keywords internal
//'
*/
// [[Rcpp::export]]
Rcpp::NumericVector pminmultinom_bonetti(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose) {
  int xlen = x.size();
  Progress prog(xlen, verbose);  // verbose = false silences the bar

  int m = prob.size();

  Rcpp::NumericVector r(xlen);
  for (int k = 0; k < xlen; k++) {
    if (Progress::check_abort()) return Rcpp::NumericVector(0);
    prog.increment();
    r(k) = equi_prob_min_leq(size, m, x(k));

    R_CheckUserInterrupt();
  }

  if (!logd) {
    return r;
  } else {
    return Rcpp::log(r);
  }
}
