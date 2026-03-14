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

//' Utility function for computing recursively a sum.
//'
//' This is an auxiliary function to compute a sum recursively.
//'
//' @param td A length-one numeric vector indicating the value to compute the
//'   function for.
//' @param n A length-one integer vector indicating the number of independent
//'   balls.
//' @param m A length-one integer vector indicating the number of independent
//'   urns/cells.
//' @param J A length-one integer vector indicating the number of largest order
//'   statistics.
//' @param sum_depth A length-one integer vector indicating the summation depth.
//' @param cur_depth A length-one integer vector indicating the current
//'   summation depth.
//' @param rangeArg Integer vector indicating the range of the summation.
//' @return A length-one numeric vector representing the summation value.
//' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
//' @seealso \code{\link{highest_order_statistics}} for computing the CDF
//'   of the sum for the first \eqn{J} largest order statistics for an
//'   equiprobable multinomial distribution.
//' @references
//'   Bonetti, M., Cirillo, P., Ogay, A. (2019), "Computing the exact
//'   distributions of some functions of the ordered multinomial counts:
//'   maximum, minimum, range and sums of order statistics", Royal Society
//'   Open Science, 6: 190198, <http://dx.doi.org/10.1098/rsos.190198>.
//' @examples
//' highest_order_statistics(6, 10, 5, 2)
//'
// [[Rcpp::export]]
double recursive_sum(const double & td, int n, int m, int J, int sum_depth,
                     int cur_depth, arma::vec rangeArg) {
  int t = (int)floor(td);

  std::vector<int> rangeArgVec(rangeArg.n_elem);
  int partial_sum = 0;
  for (arma::uword i = 0; i < rangeArg.n_elem; i++) {
    rangeArgVec[i] = (int)rangeArg(i);
    partial_sum   += rangeArgVec[i];
  }

  return recursive_sum_impl(t, n, m, J, sum_depth, cur_depth,
                            rangeArgVec, partial_sum);
}
