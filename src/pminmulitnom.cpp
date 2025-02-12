// -*- mode: C++; c-indent-level: 2; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "XOMultinom.h"

// [[Rcpp::depends("RcppArmadillo")]]

// Note: RcppExport is an alias for extern "C"

void pmin_cond(Rcpp::NumericVector indices, double x, int size, Rcpp::NumericVector prob,
  Rcpp::Environment envir) {
  Rcpp::NumericMatrix prmin = envir["prmin"];
  Rcpp::NumericVector prmin_sum = envir["prmin_sum"];
  int prmin_idx = envir["prmin_idx"];
  int indices_sum = std::accumulate(indices.begin(), indices.end(), 0.0);
  Rcpp::NumericVector indices_tmp = indices;
  indices_tmp.push_back(prmin_sum[indices_sum - (prob.size() - 2)*x]);
  prmin(prmin_idx - 1, Rcpp::_) = indices_tmp;
  prmin_idx++;
  envir["prmin"] = prmin;
  envir["prmin_idx"] = prmin_idx;
}

double px_cond_min(double x, int size, double prob) {
  double tmp = 0;

  if (size >= 0) {
    tmp = R::dbinom(x, size, prob, 0);
  }

  return tmp;
}

// Recursive function to create loops
void create_loops_min(int current_level, Rcpp::NumericVector indices, double xa, int sz,
  Rcpp::NumericVector prb, Rcpp::Environment env, int levels,
  void (*act)(Rcpp::NumericVector, double, int, Rcpp::NumericVector, Rcpp::Environment)) {
  if (current_level > levels) {
    // Base case: All loops are complete, perform the action
    act(indices, xa, sz, prb, env);
  } else {
    // Recursive case: Create the current loop and recurse
    int indices_sum = std::accumulate(indices.begin(), indices.end(), 0.0);
    Rcpp::List prx_sum = env["prx_sum"];
    Rcpp::NumericMatrix prx_mat = prx_sum[current_level - 1];
    Rcpp::NumericVector prx_tmp;
    for (int i = xa; i <= sz; i++) {
      Rcpp::List prx = env["prx"];
      if (prx[current_level - 1] != R_NilValue) {
        prx_tmp = prx[current_level - 1];
      }
      prx_tmp.push_back(prx_mat(i - xa, indices_sum - (current_level - 1)*xa));
      prx[current_level - 1] = prx_tmp;
      env["prx"] = prx;
      Rcpp::NumericVector indices_next = indices;
      indices_next.push_back(i);
      create_loops_min(current_level + 1, indices_next, xa, sz, prb, env, levels, act);
    }
  }
}

// Define a function to dynamically create and execute nested loops
void dynamic_nested_loops_min(int levels,
  void (*action)(Rcpp::NumericVector, double, int, Rcpp::NumericVector, Rcpp::Environment),
  double x, int size, Rcpp::NumericVector prob, Rcpp::Environment envir) {
  int one = 1;
  Rcpp::NumericVector null;
  create_loops_min(one, null, x, size, prob, envir, levels, action);
}

// [[Rcpp::export]]
double pminmultinom_C_one(const double& x, const int& size, const Rcpp::NumericVector& prob,
  Rcpp::Environment& this_env) {
  int k = prob.size();
  int km2 = k - 2;
  double x_tmp = x + 1.0;  // needed to convert P(min >= x) to P(min <= x)
  Rcpp::NumericVector prob_c = cumsum_rcpp(prob);

  if (x >= 0 && x < size) {
    int prmin_idx = 1;
    this_env["prmin_idx"] = prmin_idx;
    int prmin_size = std::pow(size - x, km2);
    if (prmin_size == 0) {
      prmin_size = 4;
    }
    Rcpp::NumericMatrix prmin(prmin_size, (k - 1));
    this_env["prmin"] = prmin;

    Rcpp::NumericVector prmin_sum(km2*(size - x_tmp) + 1);
    double tmp;
    for (int i = (km2*x_tmp); i <= (km2*size); i++) {
      if ((x_tmp <= (size - i)/2) && (x_tmp >= 0)) {
        tmp =
          R::pbinom(size - i - x_tmp, size - i, prob[1]/prob_c[1], 1, 0) -
          R::pbinom(x_tmp - 1, size - i, prob[1]/prob_c[1], 1, 0);
      } else if (x_tmp > (size - i)/2) {
        tmp = 0;
      } else if (x_tmp < 0) {
        tmp = 1;
      }
      prmin_sum[i - km2*x_tmp] = tmp;
    }
    this_env["prmin_sum"] = prmin_sum;
    // Rcpp::print(prmin_sum); // [[DEBUG]]

    Rcpp::List prx_sum(km2);
    for (int p = 0; p < km2; p++) {
      Rcpp::NumericMatrix prx_mat(size - x_tmp + 1, p*(size - x_tmp) + 1);
      for (int i = x_tmp; i <= size; i++) {
        for (int j = (p*x_tmp); j <= (p*size); j++) {
          prx_mat(i - x_tmp, j - p*x_tmp) = px_cond_min(i, size - j, prob[k - p - 1]/prob_c[k - p - 1]);
        }
      }
      prx_sum[p] = prx_mat;
    }
    this_env["prx_sum"] = prx_sum;
    // Rcpp::print(prx_sum); // [[DEBUG]]

    Rcpp::List prx(km2);
    this_env["prx"] = prx;
    dynamic_nested_loops_min(km2, pmin_cond, x_tmp, size, prob, this_env);

    // Rcpp::print(prx); // [[DEBUG]]
    Rcpp::NumericMatrix red(Rcpp::clone(prmin));
    int idx = km2;
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
      red_df = aggregate_sum_rcpp(red_df, by_cols_df);
      red = df2nm_rcpp(red_df);
      Rcpp::NumericMatrix red_rev = flipcols_rcpp(red(Rcpp::_, Rcpp::Range(0, red.ncol() - 2)));
      for (int j = 0; j < (red.ncol() - 2); j++) {
        Rcpp::NumericVector v = red_rev(Rcpp::_, j);
        red(Rcpp::_, j) = v;
      }
      idx--;

      R_CheckUserInterrupt();
    }

    return (1.0 - red(0, 1));
  } else if (x < 0) {
    return 0;
  } else {
    return 1;
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector pminmultinom_C(const Rcpp::NumericVector& x, const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose, Rcpp::Environment& this_env, const double& tol) {
  int xlen = x.size();
  // double s = Rcpp::sum(prob);

  // if (Rcpp::any(!Rcpp::is_finite(prob))) {
  //   Rcpp::stop("probabilities must be finite");
  // }
  // if (any_sug(prob < 0)) {
  //   Rcpp::stop("probabilities must be non-negative");
  // }
  // if (s == 0) {
  //   Rcpp::stop("probabilities must be not all 0");
  // }

  // Rcpp::NumericVector prob_new = prob/s;
  // Rcpp::LogicalVector i0 = (prob_new == 0);
  // if (any_sug(i0)) {
  //   prob_new = prob[!i0];
  //   k = prob_new.size();
  // }
  
  Rcpp::NumericVector r(xlen);
  for (int m = 0; m < xlen; m++) {
    if (verbose) Rprintf("computing P(min(X1,..., Xk) <= x) for x = %.0f...\n", x(m));
    if (xlen > 1 && Rcpp::is_true(all(diff(x) == 1)) && m > 0 && abs(1 - r(m - 1)) < tol) {
      r(m) = 1;
    } else {
      r(m) = pminmultinom_C_one(x(m), size, prob, this_env);
    }

    R_CheckUserInterrupt();
  }

  if (!logd) {
    return r;
  } else {
    return Rcpp::log(r);
  }
}
