#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

/*
//' Probability mass for the maximum multinomial count (equiprobable case)
//'
//' Computes \eqn{P(\max_k N_k = c)} for multinomial cell counts with total
//' count `n` and equal category probabilities.
//'
//' @param n Integer total count.
//' @param m Integer number of classes.
//' @param c Numeric count value at which to evaluate the probability mass.
//' @return A probability between 0 and 1.
//'
//' @keywords internal
//'
*/
// [[Rcpp::export]]
double equi_prob_max_eq(int n, int m, double c) {
  if (!std::isfinite(c) || c != std::floor(c)) return 0.0;
  const int ci = static_cast<int>(c);
  if (ci < 0 || ci > n) return 0.0;
  double hi = max_order_statistic(static_cast<double>(ci), n, m);
  double lo = (ci > 0) ? max_order_statistic(static_cast<double>(ci - 1), n, m) : 0.0;
  return std::max(0.0, hi - lo);
}

/*
//' Probability that the minimum multinomial count is at most a threshold
//'   (equiprobable case)
//'
//' Computes \eqn{P(\min_k N_k \le c)} for multinomial cell counts with total
//' count `n` and equal category probabilities.
//'
//' @param n Integer total count.
//' @param m Integer number of classes.
//' @param c Numeric threshold for the minimum cell count.
//' @return A probability between 0 and 1.
//'
//' @keywords internal
//'
*/
// [[Rcpp::export]]
double equi_prob_min_leq(int n, int m, double c) {
  if (std::isnan(c)) return 0.0;
  if (c < 0.0)       return 0.0;
  if (c >= static_cast<double>(n)) return 1.0;
  return std::max(0.0, 1.0 - smallest_order_value(std::floor(c) + 1.0, n, m));
}

/*
//' Probability mass for the minimum multinomial count (equiprobable case)
//'
//' Computes \eqn{P(\min_k N_k = c)} for multinomial cell counts with total
//' count `n` and equal category probabilities.
//'
//' @param n Integer total count.
//' @param m Integer number of classes.
//' @param c Numeric count value at which to evaluate the probability mass.
//' @return A probability between 0 and 1.
//'
//' @keywords internal
//'
*/
// [[Rcpp::export]]
double equi_prob_min_eq(int n, int m, double c) {
  if (!std::isfinite(c) || c != std::floor(c)) return 0.0;
  const int ci = static_cast<int>(c);
  if (ci < 0 || ci > n) return 0.0;
  double hi = smallest_order_value(static_cast<double>(ci), n, m);
  double lo = smallest_order_value(static_cast<double>(ci + 1), n, m);
  return std::max(0.0, hi - lo);
}

/*
//' Probability mass for the multinomial count range
//'
//' Computes \eqn{P(\max_k N_k - \min_k N_k = r)} for multinomial cell counts
//' with total count `n` and equal cell probabilities.
//'
//' @param n Integer total count.
//' @param m Integer number of classes.
//' @param r Numeric count-range value at which to evaluate the probability mass.
//' @return A probability between 0 and 1.
//'
//' @keywords internal
//'
*/
// [[Rcpp::export]]
double equi_prob_range_eq(int n, int m, double r) {
  if (!std::isfinite(r) || r != std::floor(r)) return 0.0;
  const int ri = static_cast<int>(r);
  if (ri < 0 || ri > n) return 0.0;
  double hi = range_probability(static_cast<double>(ri), n, m);
  double lo = (ri > 0) ? range_probability(static_cast<double>(ri - 1), n, m) : 0.0;
  return std::max(0.0, hi - lo);
}

/*
//' Probability mass for the sum of the J largest multinomial order statistics
//'
//' Computes \eqn{P(\sum_k N_{<k>} = x)} for the sum of the J largest order
//' statistics of multinomial cell counts with total count `n` and equal
//' cell probabilities.
//'
//' @param n Integer total count.
//' @param m Integer number of classes.
//' @param r Numeric count-range value at which to evaluate the probability mass.
//' @return A probability between 0 and 1.
//'
//' @keywords internal
//'
*/
// [[Rcpp::export]]
double equi_prob_highest_eq(int n, int m, double r, int J) {
  if (!std::isfinite(r) || r != std::floor(r)) return 0.0;
  const int ri = static_cast<int>(r);
  if (ri < 0 || ri > n) return 0.0;
  double hi = highest_order_statistics(static_cast<double>(ri), n, m, J);
  double lo = (ri > 0) ? highest_order_statistics(static_cast<double>(ri - 1), n, m, J) : 0.0;
  return std::max(0.0, hi - lo);
}
