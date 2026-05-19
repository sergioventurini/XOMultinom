#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

/*
//' Range probability for multinomial distribution
//'
//' Computes the probability that the range of a multinomial sample
//' (i.e., max(X) - min(X)) is equal to a given value.
//'
//' @param x Numeric vector of range values at which to evaluate the probability.
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
Rcpp::NumericVector drangemultinom_corrado(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose) {
  int xlen = x.size();
  Progress prog(xlen, verbose);  // verbose = false silences the bar

  std::vector<double> pi = Rcpp::as<std::vector<double>>(prob);

  Rcpp::NumericVector r(xlen);
  for (int k = 0; k < xlen; k++) {
    if (Progress::check_abort()) return Rcpp::NumericVector(0);
    prog.increment();
    r(k) = prob_range_eq(size, pi, x(k));

    R_CheckUserInterrupt();
  }

  if (!logd) {
    return r;
  } else {
    return Rcpp::log(r);
  }
}
