#ifndef XOMultinom_H
#define XOMultinom_H

// #define ARMA_NO_DEBUG // this macro disables bounds checks in the Armadillo
                         // library making the code faster (but also more
                         // frail!); the suggestion is to disable it only for
                         // the final release (see
                         // http://arma.sourceforge.net/docs.html)

#include <R.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>
#include <climits>
#include <string>
#include <algorithm>

static const double machine_eps = 2.220446049250313080847e-16;
static const double log_pi = std::log(M_PI);
static const double log_2pi = std::log(2.0 * M_PI);
static const double log_two = std::log(2.0);

// MAIN FUNCTIONS -----------------------------------------------------------------------------------------------------

// DISTRIBUTION FUNCTIONS ---------------------------------------------------------------------------------------------

// MATRIX UTILITIES ---------------------------------------------------------------------------------------------------

// MCMC SIMULATION ----------------------------------------------------------------------------------------------------

// UTILITIES ----------------------------------------------------------------------------------------------------------

#endif
