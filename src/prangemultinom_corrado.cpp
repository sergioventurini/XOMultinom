#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

/*
//' CDF of the range for a multinomial distribution
//'
//' Computes the cumulative distribution function of the range
//' \eqn{R = \max(N_1, \ldots, N_k) - \min(N_1, \ldots, N_k)}
//' for a multinomial distribution with arbitrary cell probabilities.
//'
//' @param x Numeric vector of values at which to evaluate the CDF.
//' @param size Integer number of trials.
//' @param prob Numeric vector of category probabilities. Must sum to 1.
//' @param logd Logical; if \code{TRUE}, returns log-probabilities.
//' @param verbose Logical; if \code{TRUE}, displays a progress bar.
//'
//' @return Numeric vector of the same length as \code{x} containing
//'   \eqn{P(R \le x)} (or log-probabilities if \code{logd = TRUE}).
//'
//' @examples
//' prangemultinom_corrado(0:3, size = 10, prob = rep(1/5, 5),
//'                        logd = FALSE, verbose = FALSE)
//'
//' @keywords internal
//'
*/
// [[Rcpp::export]]
Rcpp::NumericVector prangemultinom_corrado(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose) {
  int xlen = x.size();
  Progress prog(xlen, verbose);  // verbose = false silences the bar

  std::vector<double> pi = Rcpp::as<std::vector<double>>(prob);

  Rcpp::NumericVector r(xlen);
  for (int k = 0; k < xlen; k++) {
    if (Progress::check_abort()) return Rcpp::NumericVector(0);
    prog.increment();
    r(k) = prob_range_leq(size, pi, x(k));

    R_CheckUserInterrupt();
  }

  if (!logd) {
    return r;
  } else {
    return Rcpp::log(r);
  }
}
