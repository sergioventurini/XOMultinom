// -*- mode: C++; c-indent-level: 2; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "XOMultinom.h"

// [[Rcpp::depends("RcppArmadillo")]]

// Note: RcppExport is an alias for extern "C"

// [[Rcpp::export]]
std::vector<std::vector<double>> computeQk_min(const double& x, const int& k, const int& size,
  const Rcpp::NumericVector& prob, const bool& verbose, const double& tol) {
  double x_tmp = x + 1.0;  // needed to convert P(min >= x) to P(min <= x)
  int x_int = static_cast<int>(x_tmp);
  int m = prob.size();
  if (k < 1 || k > m) {
    std::printf("k must be between 1 and %d.\n", static_cast<int>(m));
    return {};
  }

  std::vector<double> prob_rev(prob.begin(), prob.end());
  std::reverse(prob_rev.begin(), prob_rev.end()); // reverse the prob vector
  std::vector<double> prob_c(m);
  std::partial_sum(prob_rev.begin(), prob_rev.end(), prob_c.begin());
  std::reverse(prob_c.begin(), prob_c.end()); // reverse back to correct order

  int Qk_r = (k == 1) ? 1 : (size + 1);
  int Qk_c = (k == m) ? 1 : (size + 1);
  std::vector<std::vector<double>> Qk(Qk_r, std::vector<double>(Qk_c, 0.0));

  if (k == 1) {
    for (int s_k = 0; s_k <= size; s_k++) {
      if (s_k >= x_int) {
        Qk[0][s_k] = R::dbinom(s_k, size, prob[0], 0);
      }
    }
  }
  else if (k == m) {
    for (int s_km1 = 0; s_km1 <= size; s_km1++) {
      if (size - s_km1 >= x_int) {
        Qk[s_km1][0] = 1.0;
      }
    }
  }
  else {
    for (int s_km1 = 0; s_km1 <= size; s_km1++) {
      for (int s_k = s_km1; s_k <= size; s_k++) {
        if (s_k - s_km1 >= x_int) {
          Qk[s_km1][s_k] = R::dbinom(s_k - s_km1, size - s_km1, prob[k - 1]/prob_c[k - 1], 0);
        }
      }
    }
  }

  return Qk;
}

// [[Rcpp::export]]
double pminmultinom_corrado_one(const double& x, const int& size, const Rcpp::NumericVector& prob,
  const bool& verbose, const double& tol) {
  if (x < 0) {
    return 0.0;
  }
  else if (x >= size) {
    return 1.0;
  }
  else {
    int m = prob.size();
    std::vector<std::vector<double>> Qk_tmp = {};
    std::vector<std::vector<double>> res = computeQk_min(x, 2, size, prob, verbose, tol);
    for (int k = 3; k < m; k++) {
      Qk_tmp = computeQk_min(x, k, size, prob, verbose, tol);
      res = multiplyUpperTriangular(res, Qk_tmp);
    }
    res = multiplyMatrix(computeQk_min(x, 1, size, prob, verbose, tol), res);
    res = multiplyMatrix(res, computeQk_min(x, m, size, prob, verbose, tol));
    return (1.0 - res[0][0]);
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector pminmultinom_corrado(const Rcpp::NumericVector& x, const int& size, const Rcpp::NumericVector& prob,
  const bool& logd, const bool& verbose, const double& tol) {
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
    if (verbose) std::printf("computing P(min(X1,..., Xk) <= %.2f)...\n", x(m));
    if (xlen > 1 && Rcpp::is_true(all(diff(x) == 1)) && m > 0 && abs(1 - r(m - 1)) < tol) {
      r(m) = 1;
    } else {
      r(m) = pminmultinom_corrado_one(x(m), size, prob, verbose, tol);
    }

    R_CheckUserInterrupt();
  }

  if (!logd) {
    return r;
  } else {
    return Rcpp::log(r);
  }
}
