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

#endif
