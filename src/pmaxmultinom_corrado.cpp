#include "XOMultinom.h"

// [[Rcpp::depends("RcppArmadillo")]]

// Note: RcppExport is an alias for extern "C"

//' CDF of the maximum for a multinomial distribution
//'
//' Computes the cumulative distribution function of the maximum
//' \eqn{M = \max(X_1, \ldots, X_k)} for a multinomial distribution
//' with arbitrary cell probabilities.
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
//' @details
//' The function evaluates the probability that the maximum of multinomial
//' cell counts does not exceed the specified values. Computation is
//' performed element-wise over \code{x}.
//'
//' @examples
//' pmaxmultinom_corrado(0:3, size = 10, prob = rep(1/5, 5),
//'                      logd = FALSE, verbose = FALSE)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pmaxmultinom_corrado(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose) {
  int xlen = x.size();
  Progress prog(xlen, verbose);  // verbose = false silences the bar

  std::vector<double> pi = Rcpp::as<std::vector<double>>(prob);

  Rcpp::NumericVector r(xlen);
  for (int k = 0; k < xlen; k++) {
    // if (verbose) std::printf("computing P(max(X1,..., Xk) <= %.4g)...\n", x(k));
    if (Progress::check_abort()) return Rcpp::NumericVector(0); // Ctrl+C
    prog.increment();
    r(k) = prob_max_leq(size, pi, x(k));

    R_CheckUserInterrupt();
  }

  if (!logd) {
    return r;
  } else {
    return Rcpp::log(r);
  }
}
