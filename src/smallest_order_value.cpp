#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

//' Internal recursive summation helper
//'
//' Wrapper around the internal recursive summation implementation used by
//' \code{highest_order_statistics()}.
//'
//' @param td Numeric value used in the recursive computation. Internally, this
//'   is converted to \code{floor(td)}.
//' @param n Integer number of balls/trials.
//' @param m Integer number of urns/cells.
//' @param J Integer number of largest order statistics.
//' @param sum_depth Integer target summation depth.
//' @param cur_depth Integer current recursion depth.
//' @param rangeArg Numeric vector containing the current recursion arguments.
//'
//' @return Numeric value giving the recursive summation result.
//'
//' @keywords internal
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
