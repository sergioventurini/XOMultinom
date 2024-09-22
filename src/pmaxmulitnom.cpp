// -*- mode: C++; c-indent-level: 2; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "XOMultinom.h"

// [[Rcpp::depends("RcppArmadillo")]]

// Note: RcppExport is an alias for extern "C"

void pmax_cond(Rcpp::IntegerVector indices, double x, int size, Rcpp::NumericVector prob,
  Rcpp::Environment envir) {
  Rcpp::NumericVector prob_c = cumsum_rcpp(prob);
  double tmp;
  if ((x <= (size - Rcpp::sum(indices))) && (x >= (size - Rcpp::sum(indices))/2)) {
    tmp =
      R::pbinom(x, size - Rcpp::sum(indices), prob[1]/prob_c[1], 1, 0) -
      R::pbinom(size - Rcpp::sum(indices) - x - 1, size - Rcpp::sum(indices), prob[1]/prob_c[1], 1, 0);
  }
  else if (x > (size - Rcpp::sum(indices))) {
    tmp = 1;
  }
  else if (x < (size - Rcpp::sum(indices))/2) {
    tmp = 0;
  }
  Rcpp::NumericMatrix prmax = envir["prmax"];
  int prmax_idx = envir["prmax_idx"];
  Rcpp::NumericVector indices_tmp = Rcpp::as<Rcpp::NumericVector>(indices);
  indices_tmp.push_back(tmp);
  prmax(prmax_idx - 1, Rcpp::_) = indices_tmp;
  prmax_idx++;
  envir["prmax"] = prmax;
  envir["prmax_idx"] = prmax_idx;
}

double px_cond_max(double x, int size, double prob) {
  double tmp;

  if (size >= 0) {
    tmp = R::dbinom(x, size, prob, 0);
  } else {
    tmp = 0;
  }

  return tmp;
}

// Recursive function to create loops
void create_loops_max(int current_level, Rcpp::IntegerVector indices, double xa, int sz,
  Rcpp::NumericVector prb, Rcpp::Environment env, int levels,
  void (*act)(Rcpp::IntegerVector, double, int, Rcpp::NumericVector, Rcpp::Environment)) {
  if (current_level > levels) {
    // Base case: All loops are complete, perform the action
    act(indices, xa, sz, prb, env);
  } else {
    // Recursive case: Create the current loop and recurse
    Rcpp::NumericVector prb_c = cumsum_rcpp(prb);
    int k = prb.size();
    int idx = k - current_level + 1;
    for (int i = 0; i <= xa; i++) {
      double px = px_cond_max(i, sz - Rcpp::sum(indices), prb[idx - 1]/prb_c[idx - 1]);
      Rcpp::List prx = env["prx"];
      Rcpp::NumericVector prx_tmp;
      if (prx[current_level - 1] != R_NilValue) {
        prx_tmp = prx[current_level - 1];
      }
      prx_tmp.push_back(px);
      prx[current_level - 1] = prx_tmp;
      env["prx"] = prx;
      Rcpp::IntegerVector indices_next = indices;
      indices_next.push_back(i);
      create_loops_max(current_level + 1, indices_next, xa, sz, prb, env, levels, act);
    }
  }
}

// Define a function to dynamically create and execute nested loops
void dynamic_nested_loops_max(int levels,
  void (*action)(Rcpp::IntegerVector, double, int, Rcpp::NumericVector, Rcpp::Environment),
  double x, int size, Rcpp::NumericVector prob, Rcpp::Environment envir) {
  int one = 1;
  Rcpp::IntegerVector null;
  create_loops_max(one, null, x, size, prob, envir, levels, action);
}

// [[Rcpp::export]]
Rcpp::NumericVector pmaxmultinom_C(const Rcpp::NumericVector& x, const int& size,const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose) {
  int k = prob.size();
  int km2 = k - 2;
  int xlen = x.size();
  double s = Rcpp::sum(prob);

  if (Rcpp::any(!Rcpp::is_finite(prob))) {
    Rcpp::stop("probabilities must be finite");
  }
  if (any_sug(prob < 0)) {
    Rcpp::stop("probabilities must be non-negative");
  }
  if (s == 0) {
    Rcpp::stop("probabilities must be not all 0");
  }

  Rcpp::NumericVector prob_new = prob/s;
  Rcpp::LogicalVector i0 = (prob_new == 0);
  if (any_sug(i0)) {
    prob_new = prob[!i0];
    k = prob_new.size();
  }

  Rcpp::Environment glob_env = Rcpp::Environment::global_env();
  Rcpp::Environment this_env = glob_env.new_child(FALSE);
  Rcpp::NumericVector r(xlen);
  for (int m = 0; m < xlen; m++) {
    if (verbose) Rprintf("computing P(max(X1,..., Xk) <= x) for x = %.0f...\n", x[m]);
    if (x[m] >= 0 && x[m] < size) {
      int prmax_idx = 1;
      this_env["prmax_idx"] = prmax_idx;
      Rcpp::NumericMatrix prmax(std::pow(x[m] + 1, km2), (k - 1));
      this_env["prmax"] = prmax;
      Rcpp::List prx(km2);
      this_env["prx"] = prx;
      dynamic_nested_loops_max(km2, pmax_cond, x[m], size, prob, this_env);

      Rcpp::NumericMatrix red = prmax;
      int idx = k - 2;
      while (idx > 0) {
        // prx = this_env["prx"];
        Rcpp::NumericVector red_tmp = red(Rcpp::_, red.ncol() - 1);
        red(Rcpp::_, red.ncol() - 1) = matelmult_rcpp(red_tmp, prx[idx - 1]);
        Rcpp::List by_cols;
        if (idx > 1) {
          by_cols = list_resize_rcpp(by_cols, red.ncol() - 2);
          for (int p = 0; p < (red.ncol() - 2); p++) {
            Rcpp::NumericVector by_col_tmp = red(Rcpp::_, p);
            by_cols[p] = Rcpp::List::create(Rcpp::Named("X") = by_col_tmp);;
          }
        } else {
          by_cols = list_resize_rcpp(by_cols, 1);
          Rcpp::IntegerVector zero_col(red.nrow());
          by_cols[0] = zero_col;
        }
        Rcpp::DataFrame by_cols_df = Rcpp::as<Rcpp::DataFrame>(by_cols);
        Rcpp::DataFrame red_df = nm2df_rcpp(red(Rcpp::_, Rcpp::Range(red.ncol() - 1, red.ncol() - 1)));
        red_df = aggregate_rcpp(red_df, by_cols_df);
        red = df2nm_rcpp(red_df);
        Rcpp::NumericMatrix red_rev = flipcols_rcpp(red(Rcpp::_, Rcpp::Range(0, red.ncol() - 2)));
        for (int j = 0; j <= (red.ncol() - 2); j++) {
          Rcpp::NumericVector v = red_rev(Rcpp::_, j);
          red(Rcpp::_, j) = v;
        }
        idx--;

        R_CheckUserInterrupt();
      }
      r(m) = red(0, 1);

      R_CheckUserInterrupt();
    } else if (x[m] < 0) {
      r(m) = 0;
    } else {
      r(m) = 1;
    }
  }

  if (!logd) {
    return r;
  } else {
    return Rcpp::log(r);
  }
}
