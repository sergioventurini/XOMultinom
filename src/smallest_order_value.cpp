#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

/*
//' Survival function of the minimum for the equiprobable multinomial case
//'
//' Computes \eqn{P(\min(N_1, \ldots, N_m) \ge t)} for a multinomial random
//' vector with equal cell probabilities.
//'
//' @param td Numeric threshold value. Internally converted to \code{floor(td)}.
//' @param n Integer number of trials.
//' @param m Integer number of cells.
//'
//' @return Numeric value giving \eqn{P(\min \ge t)}.
//'
//' @keywords internal
//'
*/
// [[Rcpp::export]]
double smallest_order_value(const double & td, int n, int m) {
  int t = (int)floor(td);

  if (t > (int)floor((double)n / m)) {
    return 0.0;
  }
  if (t == 0) {
    return 1.0;
  }

  // max_for_min now memoizes (t_max, n, m, t) triples, so repeated calls
  // with the same arguments (e.g. evaluating the CDF at many t values) are
  // effectively free after the first computation.
  return max_for_min((double)n, n, m, t);
}
