// -*- mode: C++; c-indent-level: 2; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "XOMultinom.h"

// [[Rcpp::depends("RcppArmadillo")]]

// Note: RcppExport is an alias for extern "C"

void pmin_cond(std::vector<double> indices, double x,
  std::unique_ptr<double>& res, std::unique_ptr<std::vector<double>>& prmin_sum,
  std::vector<std::unique_ptr<std::vector<std::vector<double>>>>& prx_sum) {
  int km2 = indices.size();
  int indices_sum = std::accumulate(indices.begin(), indices.end(), 0);
  double tmp = (*prmin_sum)[indices_sum - km2*x];
  int ix, indices_sum_sub = 0;
  for (int j = 0; j < km2; j++) {
    ix = static_cast<int>(indices[j]);
    indices_sum_sub = (j == 0 ? 0 : std::accumulate(indices.begin(), indices.begin() + j, 0));
    tmp *= (*prx_sum[j])[ix - x][indices_sum_sub - j*x];
  }
  *res += tmp;
}

double px_cond_min(double x, int size, double prob) {
  double tmp = 0;

  if (size >= 0) {
    tmp = R::dbinom(x, size, prob, 0);
  }

  return tmp;
}

// Recursive function to create loops
void create_loops_min(int current_level, std::vector<double> indices, double xa, int sz,
  Rcpp::NumericVector prb, int levels,
  void (*act)(std::vector<double>, double, std::unique_ptr<double>&,
    std::unique_ptr<std::vector<double>>&,
    std::vector<std::unique_ptr<std::vector<std::vector<double>>>>&), Progress& prg,
  std::unique_ptr<double>& res, std::unique_ptr<std::vector<double>>& prmin_sum,
  std::vector<std::unique_ptr<std::vector<std::vector<double>>>>& prx_sum,
  const double& tol) {
  if (*res == 1.0) return;
  if (current_level > levels) {
    // base case: all loops are complete, perform the action
    if (!Progress::check_abort()) {
      prg.increment(); // update progress
      act(indices, xa, res, prmin_sum, prx_sum);
      if (*res > 1.0 - tol) {
        *res = 1.0;
        return;
      }
    }
  } else {
    // recursive case: create the current loop and recurse
    for (int i = xa; i <= sz; i++) {
      indices[current_level - 1] = i;
      create_loops_min(current_level + 1, indices, xa, sz, prb, levels, act, prg, res, prmin_sum, prx_sum, tol);

      R_CheckUserInterrupt();
    }
  }
}

// Define a function to dynamically create and execute nested loops
void dynamic_nested_loops_min(int levels,
  void (*action)(std::vector<double>, double, std::unique_ptr<double>&,
    std::unique_ptr<std::vector<double>>&,
    std::vector<std::unique_ptr<std::vector<std::vector<double>>>>&), double x, int size,
  Rcpp::NumericVector prob, bool verbose,
  std::unique_ptr<double>& res, std::unique_ptr<std::vector<double>>& prmin_sum,
  std::vector<std::unique_ptr<std::vector<std::vector<double>>>>& prx_sum,
  const double& tol) {
  int one = 1;
  std::vector<double> indices(prob.size() - 2, -999);
  ETAProgressBar pb;
  Progress prg(std::pow((size - x + 1), (prob.size() - 2)), verbose, pb); // create the progress monitor
  create_loops_min(one, indices, x, size, prob, levels, action, prg, res, prmin_sum, prx_sum, tol);
}

// [[Rcpp::export]]
double pminmultinom_C_one(const double& x, const int& size, const Rcpp::NumericVector& prob,
  const bool& verbose, const double& tol) {
  int k = prob.size();
  int km2 = k - 2;
  double x_tmp = x + 1.0;  // needed to convert P(min >= x) to P(min <= x)
  Rcpp::NumericVector prob_c = cumsum_rcpp(prob);

  auto res = std::make_unique<double>(0.0);

  if (x >= 0 && x < size) {
    auto prmin_sum = std::make_unique<std::vector<double>>(km2*(size - x_tmp) + 1);
    double tmp = 0;
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
      (*prmin_sum)[i - km2*x_tmp] = tmp;
    }

    std::vector<std::unique_ptr<std::vector<std::vector<double>>>> prx_sum(km2);
    for (int p = 0; p < km2; p++) {
      prx_sum[p] = std::make_unique<std::vector<std::vector<double>>>(size - x_tmp + 1, std::vector<double>(p*(size - x_tmp) + 1));
      for (int i = x_tmp; i <= size; i++) {
        for (int j = (p*x_tmp); j <= (p*size); j++) {
          (*prx_sum[p])[i - x_tmp][j - p*x_tmp] = px_cond_min(i, size - j, prob[k - p - 1]/prob_c[k - p - 1]);
        }
      }
    }

    dynamic_nested_loops_min(km2, pmin_cond, x_tmp, size, prob, verbose, res, prmin_sum, prx_sum, tol);

    return (1.0 - *res);
  } else if (x < 0) {
    return 0;
  } else {
    return 1;
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector pminmultinom_C(const Rcpp::NumericVector& x, const int& size, const Rcpp::NumericVector& prob,
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
      r(m) = pminmultinom_C_one(x(m), size, prob, verbose, tol);
    }

    R_CheckUserInterrupt();
  }

  if (!logd) {
    return r;
  } else {
    return Rcpp::log(r);
  }
}
