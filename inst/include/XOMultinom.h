#ifndef XOMULTINOM_H
#define XOMULTINOM_H

#include <Rcpp.h>
#include <R_ext/Utils.h>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <functional>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <progress.hpp>

static const double machine_eps = 2.220446049250313080847e-16;
static const double pi = 3.14159265358979323846;
static const double log_pi  = std::log(pi);
static const double log_2pi = std::log(2.0 * pi);
static const double log_two = std::log(2.0);

// ---------------------------------------------------------------------------
// Hash functors for memoization caches
// ---------------------------------------------------------------------------
// Hash for std::tuple<int, int, int>
struct TupleHash3 {
  std::size_t operator()(const std::tuple<int, int, int>& tpl) const noexcept {
    std::size_t h = 0;
    h ^= std::hash<int>{}(std::get<0>(tpl)) + 0x9e3779b9u + (h << 6) + (h >> 2);
    h ^= std::hash<int>{}(std::get<1>(tpl)) + 0x9e3779b9u + (h << 6) + (h >> 2);
    h ^= std::hash<int>{}(std::get<2>(tpl)) + 0x9e3779b9u + (h << 6) + (h >> 2);
    return h;
  }
};

// Hash for std::tuple<int, int, int, int>
struct TupleHash4 {
  std::size_t operator()(const std::tuple<int, int, int, int>& tpl) const noexcept {
    std::size_t h = 0;
    h ^= std::hash<int>{}(std::get<0>(tpl)) + 0x9e3779b9u + (h << 6) + (h >> 2);
    h ^= std::hash<int>{}(std::get<1>(tpl)) + 0x9e3779b9u + (h << 6) + (h >> 2);
    h ^= std::hash<int>{}(std::get<2>(tpl)) + 0x9e3779b9u + (h << 6) + (h >> 2);
    h ^= std::hash<int>{}(std::get<3>(tpl)) + 0x9e3779b9u + (h << 6) + (h >> 2);
    return h;
  }
};

// ---------------------------------------------------------------------------
// FUNCTIONS FROM BONETTI ET AL. (2019)
// ---------------------------------------------------------------------------
double max_order_statistic(const double & td, int n, int m);
double recursive_sum(const double & td, int n, int m, int J, int sum_depth,
                     int cur_depth, Rcpp::NumericVector rangeArg);
double highest_order_statistics(const double & td, int n, int m, int J);
double max_for_min(const double & t_max, int n, int m, int t);
double smallest_order_value(const double & td, int n, int m);
double max_for_range(const double & t_max, int n, int m,
  Rcpp::NumericVector prev, int t);
double range_probability(const double & td, int n, int m);

double max_for_range_impl(double t_max, int n, int m,
                          int orig_t_max, int orig_n, int orig_m, int t);
double recursive_sum_impl(int t, int n, int m, int J, int sum_depth,
                          int cur_depth, std::vector<int>& rangeArg,
                          int partial_sum);

double equi_prob_max_eq(int n, int m, double c);

double equi_prob_min_leq(int n, int m, double c);
double equi_prob_min_eq(int n, int m, double c);

double equi_prob_range_eq(int n, int m, double r);

double equi_prob_highest_eq(int n, int m, double r, int J);

Rcpp::NumericVector pmaxmultinom_bonetti(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose);
Rcpp::NumericVector pminmultinom_bonetti(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose);
Rcpp::NumericVector prangemultinom_bonetti(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose);
Rcpp::NumericVector pJlargemultinom_bonetti(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob, const int& J,
  const bool& logd, const bool& verbose);
Rcpp::NumericVector dmaxmultinom_bonetti(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose);
Rcpp::NumericVector dminmultinom_bonetti(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose);
Rcpp::NumericVector drangemultinom_bonetti(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose);
Rcpp::NumericVector dJlargemultinom_bonetti(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob, const int& J,
  const bool& logd, const bool& verbose);

// ---------------------------------------------------------------------------
// FUNCTIONS FROM CORRADO (2011)
// ---------------------------------------------------------------------------
double prob_max_leq(int n, const std::vector<double>& pi, double c);
double prob_max_gt(int n, const std::vector<double>& pi, double c);
double prob_max_eq(int n, const std::vector<double>& pi, double c);

double prob_min_geq(int n, const std::vector<double>& pi, double c);
double prob_min_lt(int n, const std::vector<double>& pi, double c);
double prob_min_leq(int n, const std::vector<double>& pi, double c);
double prob_min_eq(int n, const std::vector<double>& pi, double c);

double prob_joint(int n, const std::vector<double>& pi, double a, double b);

double prob_range_lt(int n, const std::vector<double>& pi, double r);
double prob_range_leq(int n, const std::vector<double>& pi, double r);
double prob_range_geq(int n, const std::vector<double>& pi, double r);
double prob_range_gt(int n, const std::vector<double>& pi, double r);
double prob_range_eq(int n, const std::vector<double>& pi, double r);

Rcpp::NumericVector pmf_max_range(int n, const std::vector<double>& pi,
                                  double c_lo, double c_hi);
Rcpp::NumericVector pmf_min_range(int n, const std::vector<double>& pi,
                                  double c_lo, double c_hi);
Rcpp::NumericVector pmf_range_range(int n, const std::vector<double>& pi,
                                    double r_lo, double r_hi);

Rcpp::NumericVector pmaxmultinom_corrado(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose);
Rcpp::NumericVector pminmultinom_corrado(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose);
Rcpp::NumericVector prangemultinom_corrado(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose);
Rcpp::NumericVector dmaxmultinom_corrado(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose);
Rcpp::NumericVector dminmultinom_corrado(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose);
Rcpp::NumericVector drangemultinom_corrado(const Rcpp::NumericVector& x,
  const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose);

#endif
