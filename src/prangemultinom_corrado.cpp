#include "XOMultinom.h"

// [[Rcpp::depends("RcppArmadillo")]]

// Note: RcppExport is an alias for extern "C"

// [[Rcpp::export]]
Rcpp::NumericVector prangemultinom_corrado(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose) {
  int xlen = x.size();

  std::vector<double> pi = Rcpp::as<std::vector<double>>(prob);

  Rcpp::NumericVector r(xlen);
  for (int k = 0; k < xlen; k++) {
    if (verbose) std::printf("computing P(range(X1,..., Xk) <= %.4g)...\n", x(k));
    r(k) = prob_range_leq(size, pi, x(k));

    R_CheckUserInterrupt();
  }

  if (!logd) {
    return r;
  } else {
    return Rcpp::log(r);
  }
}
