// -*- mode: C++; c-indent-level: 2; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "XOMultinom.h"

// [[Rcpp::depends("RcppArmadillo")]]

// Note: RcppExport is an alias for extern "C"

// [[Rcpp::export]]
double prangemultinom_corrado_one(const double& x, const int& size, const Rcpp::NumericVector& prob,
  const bool& verbose, const double& tol) {
  int x_tmp = static_cast<int>(x + 1);
  if (x < 0) {
    return 0.0;
  }
  else if (x >= size) {
    return 1.0;
  }
  else {
    int m = prob.size();

    std::vector<std::unique_ptr<std::vector<std::vector<double>>>> Qk(m);
    for (int k = 1; k <= m; k++) {
      Qk[k - 1] = std::make_unique<std::vector<std::vector<double>>>(size + 1, std::vector<double>(size + 1));
      *Qk[k - 1] = computeQk_full(k, size, prob, verbose, tol);
    }

    std::vector<std::vector<double>> Qk_tmp = {};
    std::vector<std::vector<double>> res_partI = {};
    double res = 0.0;
    for (int h = 0; h <= (size - x_tmp + 1); h++) {
      res_partI = computeQk_culled(*Qk[1], h + x_tmp - 1, h, 2, size, prob, verbose, tol);
      for (int k = 3; k < m; k++) {
        Qk_tmp = computeQk_culled(*Qk[k - 1], h + x_tmp - 1, h, k, size, prob, verbose, tol);
        res_partI = multiplyUpperTriangular(res_partI, Qk_tmp);
      }
      res_partI = multiplyMatrix(
        computeQk_culled(*Qk[0], h + x_tmp - 1, h, 1, size, prob, verbose, tol), res_partI);
      res_partI = multiplyMatrix(res_partI,
        computeQk_culled(*Qk[m - 1], h + x_tmp - 1, h, m, size, prob, verbose, tol));
      res += res_partI[0][0];
    }

    std::vector<std::vector<double>> res_partII = {};
    for (int h = 0; h <= (size - x_tmp); h++) {
      res_partII = computeQk_culled(*Qk[1], h + x_tmp - 1, h + 1, 2, size, prob, verbose, tol);
      for (int k = 3; k < m; k++) {
        Qk_tmp = computeQk_culled(*Qk[k - 1], h + x_tmp - 1, h + 1, k, size, prob, verbose, tol);
        res_partII = multiplyUpperTriangular(res_partII, Qk_tmp);
      }
      res_partII = multiplyMatrix(
        computeQk_culled(*Qk[0], h + x_tmp - 1, h + 1, 1, size, prob, verbose, tol), res_partII);
      res_partII = multiplyMatrix(res_partII,
        computeQk_culled(*Qk[m - 1], h + x_tmp - 1, h + 1, m, size, prob, verbose, tol));
      res -= res_partII[0][0];
    }

    return res;
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector prangemultinom_corrado(const Rcpp::NumericVector& x, const int& size,
  const Rcpp::NumericVector& prob, const bool& logd, const bool& verbose, const double& tol) {
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
    if (verbose) std::printf("computing P(range(X1,..., Xk) <= %.2f)...\n", x(m));
    if (xlen > 1 && Rcpp::is_true(all(diff(x) == 1)) && m > 0 && abs(1 - r(m - 1)) < tol) {
      r(m) = 1;
    } else {
      r(m) = prangemultinom_corrado_one(x(m), size, prob, verbose, tol);
    }

    R_CheckUserInterrupt();
  }

  if (!logd) {
    return r;
  } else {
    return Rcpp::log(r);
  }
}
