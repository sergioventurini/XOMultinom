// -*- mode: C++; c-indent-level: 2; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "XOMultinom.h"

// [[Rcpp::depends("RcppArmadillo")]]

// Note: RcppExport is an alias for extern "C"

void pmax_cond(std::vector<double> indices,
  std::unique_ptr<double>& res, std::unique_ptr<std::vector<double>>& prmax_sum,
  std::vector<std::unique_ptr<std::vector<std::vector<double>>>>& prx_sum) {
  int km2 = indices.size();
  int indices_sum = std::accumulate(indices.begin(), indices.end(), 0);
  double tmp = (*prmax_sum)[indices_sum];
  int ix, indices_sum_sub = 0;
  for (int j = 0; j < km2; j++) {
    ix = static_cast<int>(indices[j]);
    indices_sum_sub = (j == 0 ? 0 : std::accumulate(indices.begin(), indices.begin() + j, 0));
    if (tmp > 0.0) {
      tmp *= (*prx_sum[j])[ix][indices_sum_sub];
    }
    else {
      tmp = 0.0;
      break;
    }
  }
  *res += tmp;
}

double px_cond_max(double x, int size, double prob) {
  double tmp = (size >= 0) ? R::dbinom(x, size, prob, 0) : 0;

  // if (size >= 0) {
  //   tmp = R::dbinom(x, size, prob, 0);
  // }

  return tmp;
}

// Recursive function to create loops
void create_loops_max(int current_level, std::vector<double> indices, double xa, int sz,
  Rcpp::NumericVector prb, int levels,
  void (*act)(std::vector<double>, std::unique_ptr<double>&,
    std::unique_ptr<std::vector<double>>&,
    std::vector<std::unique_ptr<std::vector<std::vector<double>>>>&), Progress& prg,
  std::unique_ptr<double>& res, std::unique_ptr<std::vector<double>>& prmax_sum,
  std::vector<std::unique_ptr<std::vector<std::vector<double>>>>& prx_sum,
  const double& tol) {
  if (*res == 1.0) return;
  if (current_level > levels) {
    // base case: all loops are complete, perform the action
    if (!Progress::check_abort()) {
      prg.increment(); // update progress
      act(indices, res, prmax_sum, prx_sum);
      if (*res > 1.0 - tol) {
        *res = 1.0;
        return;
      }
    }
  } else {
    // recursive case: create the current loop and recurse
    for (int i = 0; i <= xa; i++) {
      indices[current_level - 1] = i;
      create_loops_max(current_level + 1, indices, xa, sz, prb, levels, act, prg, res, prmax_sum, prx_sum, tol);

      R_CheckUserInterrupt();
    }
  }
}

// Define a function to dynamically create and execute nested loops
void dynamic_nested_loops_max(int levels,
  void (*action)(std::vector<double>, std::unique_ptr<double>&,
    std::unique_ptr<std::vector<double>>&,
    std::vector<std::unique_ptr<std::vector<std::vector<double>>>>&), double x, int size,
  Rcpp::NumericVector prob, bool verbose,
  std::unique_ptr<double>& res, std::unique_ptr<std::vector<double>>& prmax_sum,
  std::vector<std::unique_ptr<std::vector<std::vector<double>>>>& prx_sum,
  const double& tol) {
  int one = 1;
  std::vector<double> indices(prob.size() - 2, -999);
  ETAProgressBar pb;
  Progress prg(std::pow((x + 1), (prob.size() - 2)), verbose, pb); // create the progress monitor
  create_loops_max(one, indices, x, size, prob, levels, action, prg, res, prmax_sum, prx_sum, tol);
}

// [[Rcpp::export]]
double pmaxmultinom_C_one(const double& x, const int& size, const Rcpp::NumericVector& prob,
  const bool& verbose, const double& tol) {
  int k = prob.size();
  int km2 = k - 2;
  Rcpp::NumericVector prob_c = cumsum_rcpp(prob);

  auto res = std::make_unique<double>(0.0);

  if (x >= 0 && x < size) {
    auto prmax_sum = std::make_unique<std::vector<double>>(km2*x + 1, 0.0);
    double tmp = 0;
    for (int i = 0; i <= km2*x; i++) {
      if ((x <= (size - i)) && (x >= (size - i)/2)) {
        tmp =
          R::pbinom(x, size - i, prob[1]/prob_c[1], 1, 0) -
          R::pbinom(size - i - x - 1, size - i, prob[1]/prob_c[1], 1, 0);
      } else if (x > (size - i)) {
        tmp = 1;
      } else if (x < (size - i)/2) {
        tmp = 0;
      }
      (*prmax_sum)[i] = tmp;
    }

    std::vector<std::unique_ptr<std::vector<std::vector<double>>>> prx_sum(km2);
    for (int p = 0; p < km2; p++) {
      prx_sum[p] = std::make_unique<std::vector<std::vector<double>>>(x + 1, std::vector<double>(p * x + 1));
      for (int i = 0; i <= x; i++) {
        for (int j = 0; j < (p*x + 1); j++) {
          (*prx_sum[p])[i][j] = px_cond_max(i, size - j, prob[k - p - 1]/prob_c[k - p - 1]);
        }
      }
    }

    dynamic_nested_loops_max(km2, pmax_cond, x, size, prob, verbose, res, prmax_sum, prx_sum, tol);

    return *res;
  } else if (x < 0) {
    return 0.0;
  } else {
    return 1.0;
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector pmaxmultinom_C(const Rcpp::NumericVector& x, const int& size, const Rcpp::NumericVector& prob,
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
    if (verbose) std::printf("computing P(max(X1,..., Xk) <= %.2f)...\n", x(m));
    if (xlen > 1 && Rcpp::is_true(all(diff(x) == 1)) && m > 0 && abs(1 - r(m - 1)) < tol) {
      r(m) = 1;
    } else {
      r(m) = pmaxmultinom_C_one(x(m), size, prob, verbose, tol);
    }

    R_CheckUserInterrupt();
  }

  if (!logd) {
    return r;
  } else {
    return Rcpp::log(r);
  }
}
