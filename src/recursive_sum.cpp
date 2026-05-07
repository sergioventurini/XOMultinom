#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

double recursive_sum_impl(int t, int n, int m, int J, int sum_depth,
                          int cur_depth, std::vector<int>& rangeArg,
                          int partial_sum) {
  double S = 0.0;

  if (cur_depth <= sum_depth) {
    int lo, hi;
    if (cur_depth == 1) {
      lo = (int)floor((double)t / J) + 1;
      hi = t - sum_depth + 1;
    } else {
      lo = (int)floor((double)(t - partial_sum) / (J - cur_depth + 1)) + 1;
      hi = std::min(rangeArg[cur_depth - 2], t - partial_sum);
    }

    for (int r = lo; r <= hi; r++) {
      rangeArg.push_back(r);
      S += recursive_sum_impl(t, n, m, J, sum_depth,
                              cur_depth + 1, rangeArg,
                              partial_sum + r);
      rangeArg.pop_back();
      R_CheckUserInterrupt();
    }

  } else {
    const int n_rem = n - partial_sum;   // n - sum(rangeArg)

    const double prob_arg = std::floor((double)(t - partial_sum) / (J - sum_depth));
    const double temp_p   = max_order_statistic(prob_arg, n_rem, m - sum_depth);

    const double common_term = lgamma(n + 1) + lgamma(m + 1) - (double)n * std::log((double)m);

    const double log_rem = (n_rem > 0 && m != sum_depth)
                             ? (double)n_rem * std::log((double)(m - sum_depth))
                             : 0.0;

    double coef = log_rem - lgamma(m - sum_depth + 1) - lgamma(n_rem + 1);

    for (int k = 0; k < (int)rangeArg.size(); k++) {
      coef -= lgamma((double)rangeArg[k] + 1);
    }

    std::unordered_map<int, int> freq;
    freq.reserve(rangeArg.size());
    for (int val : rangeArg) freq[val]++;
    for (auto& kv : freq) {
      coef -= lgamma((double)kv.second + 1);
    }

    S = temp_p * std::exp(common_term + coef);
  }

  return S;
}

/*
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
//'
*/
// [[Rcpp::export]]
double recursive_sum(const double & td, int n, int m, int J, int sum_depth,
                     int cur_depth, Rcpp::NumericVector rangeArg) {
  int t = (int)floor(td);

  std::vector<int> rangeArgVec(rangeArg.size());
  int partial_sum = 0;
  for (R_xlen_t i = 0; i < rangeArg.size(); i++) {
    rangeArgVec[i] = (int)rangeArg[i];
    partial_sum   += rangeArgVec[i];
  }

  return recursive_sum_impl(t, n, m, J, sum_depth, cur_depth,
                            rangeArgVec, partial_sum);
}
