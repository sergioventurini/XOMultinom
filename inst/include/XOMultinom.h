#ifndef XOMULTINOM_H
#define XOMULTINOM_H

// #define ARMA_NO_DEBUG // this macro disables bounds checks in the Armadillo
                         // library making the code faster (but also more
                         // frail!); the suggestion is to disable it only for
                         // the final release (see
                         // http://arma.sourceforge.net/docs.html)

#include <R.h>
#include <R_ext/Utils.h>
#include <RcppArmadillo.h>
#include <cmath>
#include <climits>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "progress.hpp"
#include "progbar.h"

static const double machine_eps = 2.220446049250313080847e-16;
static const double log_pi = std::log(M_PI);
static const double log_2pi = std::log(2.0 * M_PI);
static const double log_two = std::log(2.0);

double max_order_statistic_C(const double & td, int n, int m);
double recursive_sum_C(const double & td, int n, int m, int J, int sum_depth, int cur_depth, arma::vec rangeArg);
double highest_order_statistics_C(const double & td, int n, int m, int J);
double max_for_min_C(const double & t_max, int n, int m, int t);
double smallest_order_value_C(const double & td, int n, int m);
double max_for_range_C(const double & t_max, int n, int m, arma::vec prev, int t);
double range_probability_C(const double & td, int n, int m);

void pmax_cond(std::vector<double> indices, std::unique_ptr<double>& res,
  std::unique_ptr<std::vector<double>>& prmax_sum,
  std::vector<std::unique_ptr<std::vector<std::vector<double>>>>& prx_sum);
double px_cond(double x, int size, double prob);
void create_loops_max(int current_level, std::vector<double> indices, double xa, int sz,
  Rcpp::NumericVector prb, int levels,
  double (*act)(std::vector<double>, std::unique_ptr<double>&,
    std::unique_ptr<std::vector<double>>&,
    std::vector<std::unique_ptr<std::vector<std::vector<double>>>>&),
  Progress& prg, std::unique_ptr<double>& res, std::unique_ptr<std::vector<double>>& prmax_sum,
  std::vector<std::unique_ptr<std::vector<std::vector<double>>>>& prx_sum,
  const double& tol);
void dynamic_nested_loops_max(int levels,
  double (*action)(std::vector<double>, std::unique_ptr<double>&,
    std::unique_ptr<std::vector<double>>&,
    std::vector<std::unique_ptr<std::vector<std::vector<double>>>>&), double x, int size,
  Rcpp::NumericVector prob, bool verbose, std::unique_ptr<double>& res,
  std::unique_ptr<std::vector<double>>& prmax_sum,
  std::vector<std::unique_ptr<std::vector<std::vector<double>>>>& prx_sum,
  const double& tol);
double pmaxmultinom_C_one(const double& x, const int& size, const Rcpp::NumericVector& prob,
  const bool& verbose, const double& tol);
Rcpp::NumericVector pmaxmultinom_C(const Rcpp::NumericVector& x, const int& size, const Rcpp::NumericVector& probs, const bool& logd, const bool& verbose, const double& tol);

void pmin_cond(std::vector<double> indices, double x, std::unique_ptr<double>& res,
  std::unique_ptr<std::vector<double>>& prmin_sum,
  std::vector<std::unique_ptr<std::vector<std::vector<double>>>>& prx_sum);
double px_cond_min(double x, int size, double prob);
void create_loops_min(int current_level, std::vector<double> indices, double xa, int sz,
  Rcpp::NumericVector prb, int levels,
  double (*act)(std::vector<double>, double, std::unique_ptr<double>&,
    std::unique_ptr<std::vector<double>>&,
    std::vector<std::unique_ptr<std::vector<std::vector<double>>>>&),
  Progress& prg, std::unique_ptr<double>& res, std::unique_ptr<std::vector<double>>& prmin_sum,
  std::vector<std::unique_ptr<std::vector<std::vector<double>>>>& prx_sum,
  const double& tol);
void dynamic_nested_loops_min(int levels,
  double (*action)(std::vector<double>, double, std::unique_ptr<double>&,
    std::unique_ptr<std::vector<double>>&,
    std::vector<std::unique_ptr<std::vector<std::vector<double>>>>&), double x, int size,
  Rcpp::NumericVector prob, bool verbose, std::unique_ptr<double>& res,
  std::unique_ptr<std::vector<double>>& prmin_sum,
  std::vector<std::unique_ptr<std::vector<std::vector<double>>>>& prx_sum,
  const double& tol);
double pminmultinom_C_one(const double& x, const int& size, const Rcpp::NumericVector& prob,
  const bool& verbose, const double& tol);
Rcpp::NumericVector pminmultinom_C(const Rcpp::NumericVector& x, const int& size, const Rcpp::NumericVector& probs, const bool& logd, const bool& verbose, const double& tol);

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

#endif
