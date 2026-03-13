#ifndef XOMULTINOM_H
#define XOMULTINOM_H

// #define ARMA_NO_DEBUG // this macro disables bounds checks in the Armadillo
                         // library making the code faster (but also more
                         // frail!); the suggestion is to disable it only for
                         // the final release (see
                         // http://arma.sourceforge.net/docs.html)

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <R.h>
#include <R_ext/Utils.h>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include "progress.hpp"
#include "progbar.h"

static const double machine_eps = 2.220446049250313080847e-16;
static const double log_pi = std::log(M_PI);
static const double log_2pi = std::log(2.0 * M_PI);
static const double log_two = std::log(2.0);

// FUNCTIONS FROM BONETTI ET AL. (2019) ---------------------------------------
double max_order_statistic(const double & td, int n, int m);
double recursive_sum(const double & td, int n, int m, int J, int sum_depth, int cur_depth, arma::vec rangeArg);
double highest_order_statistics(const double & td, int n, int m, int J);
double max_for_min(const double & t_max, int n, int m, int t);
double smallest_order_value(const double & td, int n, int m);
double max_for_range(const double & t_max, int n, int m, arma::vec prev, int t);
double range_probability(const double & td, int n, int m);

// FUNCTIONS FROM CORRADO (2011) -----------------------------------------------
double prob_max_leq(int n, const std::vector<double>& pi, double c);
double prob_min_geq(int n, const std::vector<double>& pi, double c);
double prob_joint(int n, const std::vector<double>& pi, double a, double b);
double prob_max_gt(int n, const std::vector<double>& pi, double c);
double prob_min_lt(int n, const std::vector<double>& pi, double c);
double prob_min_leq(int n, const std::vector<double>& pi, double c);
double prob_max_eq(int n, const std::vector<double>& pi, double c);
double prob_min_eq(int n, const std::vector<double>& pi, double c);
Rcpp::NumericVector pmf_max_range(int n, const std::vector<double>& pi,
                                  double c_lo, double c_hi);
Rcpp::NumericVector pmf_min_range(int n, const std::vector<double>& pi,
                                  double c_lo, double c_hi);
double prob_range_lt(int n, const std::vector<double>& pi, double r);
double prob_range_leq(int n, const std::vector<double>& pi, double r);
double prob_range_geq(int n, const std::vector<double>& pi, double r);
double prob_range_gt(int n, const std::vector<double>& pi, double r);
double prob_range_eq(int n, const std::vector<double>& pi, double r);
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

// UTILITY FUNCTIONS ----------------------------------------------------------
bool any_sug(Rcpp::LogicalVector x);
Rcpp::NumericVector cumsum_rcpp(Rcpp::NumericVector x);
Rcpp::NumericVector matelmult_rcpp(Rcpp::NumericVector v1, Rcpp::NumericVector v2);
Rcpp::NumericVector rev_rcpp(Rcpp::NumericVector x);
Rcpp::DataFrame aggregate_sum_rcpp(Rcpp::DataFrame x, Rcpp::List by);
Rcpp::List list_resize_rcpp(const Rcpp::List& x, int newsize);
Rcpp::DataFrame nm2df_rcpp(Rcpp::NumericMatrix x);
Rcpp::NumericMatrix df2nm_rcpp(Rcpp::DataFrame x);
Rcpp::NumericMatrix flipcols_rcpp(Rcpp::NumericMatrix x);
void print_umap(std::unordered_map<int, double> myMap);
std::vector<std::vector<double>> addMatrix(const std::vector<std::vector<double>>& A,
  const std::vector<std::vector<double>>& B);
std::vector<std::vector<double>> multiplyMatrix(const std::vector<std::vector<double>>& A,
  const std::vector<std::vector<double>>& B);
std::vector<std::vector<double>> multiplyUpperTriangular(const std::vector<std::vector<double>>& A,
  const std::vector<std::vector<double>>& B);
Rcpp::NumericMatrix vector2D_2_NM(std::vector<std::vector<double>> mat);
arma::mat vec_2_armaMat(std::vector<std::vector<double>> x);
std::vector<std::vector<double>> armaMat_2_vec(const arma::mat& mat);
void printVector(const std::unique_ptr<std::vector<double>>& vecPtr);
void printComplexVector(const std::vector<std::unique_ptr<std::vector<std::vector<double>>>>& complexVec);
void printMatrix(const std::vector<std::vector<double>>& matrix);

#endif
