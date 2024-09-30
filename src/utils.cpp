// -*- mode: C++; c-indent-level: 2; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <R_ext/Utils.h>

// [[Rcpp::depends("RcppArmadillo")]]

// Note: RcppExport is an alias for extern "C"

// [[Rcpp::export]]
bool any_sug(Rcpp::LogicalVector x) {
   // Note the use of is_true to return a bool type
   return Rcpp::is_true(any(x == TRUE));
}

// [[Rcpp::export]]
Rcpp::NumericVector cumsum_rcpp(Rcpp::NumericVector x) {
  // initialize the result vector
  Rcpp::NumericVector res(x.size());
  std::partial_sum(x.begin(), x.end(), res.begin());
  return res;
}

// [[Rcpp::export]]
Rcpp::NumericVector matelmult_rcpp(Rcpp::NumericVector v1, Rcpp::NumericVector v2) {
  int n = v1.size();
  Rcpp::NumericVector res(n);
  for (int i = 0; i < n; i++) {
    res(i) = v1(i) * v2(i);
  }
  return res;
}

// [[Rcpp::export]]
Rcpp::NumericVector rev_rcpp(Rcpp::NumericVector x) {
  Rcpp::NumericVector rev_x = Rcpp::clone<Rcpp::NumericVector>(x);
  std::reverse(rev_x.begin(), rev_x.end());
  return rev_x;
}

// [[Rcpp::export]]
Rcpp::DataFrame aggregate_rcpp(Rcpp::DataFrame x, Rcpp::List by) {
  Rcpp::Function f("aggregate.data.frame");

  return f(Rcpp::Named("x") = x, Rcpp::Named("by") = by, Rcpp::Named("FUN") = "sum");
}

// [[Rcpp::export]]
Rcpp::List list_resize_rcpp(const Rcpp::List& x, int newsize) {
  int oldsize = x.size();
  Rcpp::List y(newsize);

  for (int i = 0; i < oldsize; i++) {
    y[i] = x[i];
  }

  return y;
}

//[[Rcpp::export]]
Rcpp::DataFrame nm2df_rcpp(Rcpp::NumericMatrix x) {
  Rcpp::List y_tmp(x.ncol());
  
  for (int j = 0; j < x.ncol(); j++) {
    Rcpp::NumericVector v = x(Rcpp::_, j);
    // y_tmp[j] = v;
    y_tmp[j] = Rcpp::List::create(Rcpp::Named("V") = v);
  }
  Rcpp::DataFrame y = Rcpp::as<Rcpp::DataFrame>(y_tmp);

  return y;
}

//[[Rcpp::export]]
Rcpp::NumericMatrix df2nm_rcpp(Rcpp::DataFrame x) {
  Rcpp::NumericVector v = x[0];
  int n = v.size();
  int p = x.size();
  Rcpp::NumericMatrix y(n, p);
  
  for (int j = 0; j < p; j++) {
    Rcpp::NumericVector v = x[j];
    y(Rcpp::_, j) = v;
  }

  return y;
}

//[[Rcpp::export]]
Rcpp::NumericMatrix flipcols_rcpp(Rcpp::NumericMatrix x) {
  int n = x.nrow();
  int p = x.ncol();
  Rcpp::NumericMatrix y(n, p);
  
  for (int j = 0; j < p; j++) {
    Rcpp::NumericVector v = x(Rcpp::_, j);
    y(Rcpp::_, p - j - 1) = v;
  }

  return y;
}
