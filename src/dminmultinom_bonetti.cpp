#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

/*
//' Minimum probability for multinomial distribution (equiprobable case)
//'
//' Computes the probability mass function for the minimum component of a
//' multinomial random vector with total count `size` and category probabilities
//' `prob`, evaluated at each value in `x`.
//'
//' @param x Numeric vector of values at which to evaluate the probability mass
//'   function.
//' @param size Integer, the total number of trials.
//' @param prob Numeric vector of category probabilities.
//' @param logd Logical, if TRUE returns log-probabilities.
//' @param verbose Logical, if TRUE shows a progress bar.
//'
//' @return Numeric vector of probabilities (or log-probabilities if \code{logd = TRUE}).
//'
//' @keywords internal
//'
*/
// [[Rcpp::export]]
Rcpp::NumericVector dminmultinom_bonetti(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose) {
  int xlen = x.size();
  Progress prog(xlen, verbose);  // verbose = false silences the bar

  int m = prob.size();

  Rcpp::NumericVector r(xlen);
  for (int k = 0; k < xlen; k++) {
    if (Progress::check_abort()) return Rcpp::NumericVector(0);
    prog.increment();
    r(k) = equi_prob_min_eq(size, m, x(k));

    R_CheckUserInterrupt();
  }

  if (!logd) {
    return r;
  } else {
    return Rcpp::log(r);
  }
}
