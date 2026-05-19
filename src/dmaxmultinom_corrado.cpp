#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

/*
//' Probability mass function of the maximum of a multinomial random vector
//'
//' Computes the probability mass function for the maximum component of a
//' multinomial random vector with total count `size` and category probabilities
//' `prob`, evaluated at each value in `x`.
//'
//' @param x Numeric vector of values at which to evaluate the probability mass
//'   function.
//' @param size Integer total count of the multinomial distribution.
//' @param prob Numeric vector of category probabilities.
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
Rcpp::NumericVector dmaxmultinom_corrado(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose) {
  int xlen = x.size();
  Progress prog(xlen, verbose);  // verbose = false silences the bar

  std::vector<double> pi = Rcpp::as<std::vector<double>>(prob);

  Rcpp::NumericVector r(xlen);
  for (int k = 0; k < xlen; k++) {
    if (Progress::check_abort()) return Rcpp::NumericVector(0);
    prog.increment();
    r(k) = prob_max_eq(size, pi, x(k));

    R_CheckUserInterrupt();
  }

  if (!logd) {
    return r;
  } else {
    return Rcpp::log(r);
  }
}
