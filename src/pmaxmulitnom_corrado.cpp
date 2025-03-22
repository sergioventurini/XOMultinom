// -*- mode: C++; c-indent-level: 2; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "XOMultinom.h"

// [[Rcpp::depends("RcppArmadillo")]]

// Note: Rcpp::Export is an alias for extern "C"

double pmaxmultinom_corrado_one(
  const std::vector<std::unique_ptr<std::vector<std::vector<double>>>>& Qk,
  const double& x, const int& size, const Rcpp::NumericVector& prob,
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
    std::vector<std::vector<double>> res = computeQk_culled(*Qk[1], x, 0, 2, size, prob, verbose, tol);
    // printMatrix(res);
    for (int k = 3; k < m; k++) {
      Qk_tmp = computeQk_culled(*Qk[k - 1], x, 0, k, size, prob, verbose, tol);
      // printMatrix(Qk_tmp);
      res = multiplyUpperTriangular(res, Qk_tmp);
    }
    // printMatrix(computeQk_culled(*Qk[0], x, 0, 1, size, prob, verbose, tol));
    res = multiplyMatrix(computeQk_culled(*Qk[0], x, 0, 1, size, prob, verbose, tol), res);
    // printMatrix(computeQk_culled(*Qk[m - 1], x, 0, m, size, prob, verbose, tol));
    // printMatrix(
    //   multiplyMatrix(
    //     computeQk_culled(*Qk[m - 2], x, 0, m - 1, size, prob, verbose, tol),
    //     computeQk_culled(*Qk[m - 1], x, 0, m, size, prob, verbose, tol)
    //   )
    // );
    res = multiplyMatrix(res, computeQk_culled(*Qk[m - 1], x, 0, m, size, prob, verbose, tol));
    return res[0][0];
  }
}

// [[Rcpp::export]]
double pmaxmultinom_corrado_one_parallel(const double& x, const int& size, const Rcpp::NumericVector& prob,
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
    std::vector<std::vector<double>> res = computeQk(x, 0, 2, size, prob, verbose, tol);
    for (int k = 3; k < m; k++) {
      Qk_tmp = computeQk(x, 0, k, size, prob, verbose, tol);
      res = multiplyUpperTriangular(res, Qk_tmp);
    }
    res = multiplyMatrix(computeQk(x, 0, 1, size, prob, verbose, tol), res);
    res = multiplyMatrix(res, computeQk(x, 0, m, size, prob, verbose, tol));
    return res[0][0];
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector pmaxmultinom_corrado(const Rcpp::NumericVector& x, const int& size, const Rcpp::NumericVector& prob,
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

  int m = prob.size();

  std::vector<std::unique_ptr<std::vector<std::vector<double>>>> Qk(m);
  for (int k = 1; k <= m; k++) {
    Qk[k - 1] = std::make_unique<std::vector<std::vector<double>>>(size + 1, std::vector<double>(size + 1, 0.0));
    *Qk[k - 1] = computeQk_full(k, size, prob, verbose, tol);
  }

  Rcpp::NumericVector r(xlen);
  for (int k = 0; k < xlen; k++) {
    if (verbose) std::printf("computing P(max(X1,..., Xk) <= %.2f)...\n", x(k));
    if (xlen > 1 && Rcpp::is_true(all(diff(x) == 1)) && k > 0 && abs(1 - r(k - 1)) < tol) {
      r(k) = 1;
    } else {
      r(k) = pmaxmultinom_corrado_one(Qk, x(k), size, prob, verbose, tol);
    }

    R_CheckUserInterrupt();
  }

  if (!logd) {
    return r;
  } else {
    return Rcpp::log(r);
  }
}
