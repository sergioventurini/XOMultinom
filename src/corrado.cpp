// -*- mode: C++; c-indent-level: 2; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "XOMultinom.h"

// [[Rcpp::depends("RcppArmadillo")]]

// Note: RcppExport is an alias for extern "C"

// [[Rcpp::export]]
std::vector<std::vector<double>> computeQk_full(const int& k, const int& size,
  const Rcpp::NumericVector& prob, const bool& verbose, const double& tol) {
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
      Qk[0][s_k] = R::dbinom(s_k, size, prob[0], 0);
    }
  }
  else if (k == m) {
    for (int s_km1 = 0; s_km1 <= size; s_km1++) {
      Qk[s_km1][0] = 1.0;
    }
  }
  else {
    for (int s_km1 = 0; s_km1 <= size; s_km1++) {
      for (int s_k = s_km1; s_k <= size; s_k++) {
        Qk[s_km1][s_k] = R::dbinom(s_k - s_km1, size - s_km1, prob[k - 1]/prob_c[k - 1], 0);
      }
    }
  }

  return Qk;
}

// [[Rcpp::export]]
std::vector<std::vector<double>> computeQk_culled(std::vector<std::vector<double>> Qk,
  const double& a, const double& b, const int& k, const int& size,
  const Rcpp::NumericVector& prob, const bool& verbose, const double& tol) {
  int a_int = static_cast<int>(a);
  int b_int = static_cast<int>(b);
  int m = prob.size();

  std::vector<std::vector<double>> Qk_culled(Qk);

  if (k == 1) {
    for (int s_k = 0; s_k <= size; s_k++) {
      if (s_k > a_int || s_k < b_int) {
        Qk_culled[0][s_k] = 0.0;
      }
    }
  }
  else if (k == m) {
    for (int s_km1 = 0; s_km1 <= size; s_km1++) {
      if (size - s_km1 > a_int || size - s_km1 < b_int) {
        Qk_culled[s_km1][0] = 0.0;
      }
    }
  }
  else {
    for (int s_km1 = 0; s_km1 <= size; s_km1++) {
      for (int s_k = s_km1; s_k <= size; s_k++) {
        if (s_k - s_km1 > a_int || s_k - s_km1 < b_int) {
          Qk_culled[s_km1][s_k] = 0.0;
        }
      }
    }
  }

  return Qk_culled;
}

// [[Rcpp::export]]
std::vector<std::vector<double>> computeQk(const double& a, const double& b, const int& k,
  const int& size, const Rcpp::NumericVector& prob, const bool& verbose, const double& tol) {
  int a_int = static_cast<int>(a);
  int b_int = static_cast<int>(b);
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
      if (s_k <= a_int && s_k >= b_int) {
        Qk[0][s_k] = R::dbinom(s_k, size, prob[0], 0);
      }
    }
  }
  else if (k == m) {
    for (int s_km1 = 0; s_km1 <= size; s_km1++) {
      if (size - s_km1 <= a_int && size - s_km1 >= b_int) {
        Qk[s_km1][0] = 1.0;
      }
    }
  }
  else {
    for (int s_km1 = 0; s_km1 <= size; s_km1++) {
      for (int s_k = s_km1; s_k <= size; s_k++) {
        if (s_k - s_km1 <= a_int && s_k - s_km1 >= b_int) {
          Qk[s_km1][s_k] = R::dbinom(s_k - s_km1, size - s_km1, prob[k - 1]/prob_c[k - 1], 0);
        }
      }
    }
  }

  return Qk;
}

// The following functions are not needed any more (to delete in future)
// [[Rcpp::export]]
std::vector<std::vector<double>> computeQk_max(const double& x, const int& k, const int& size,
  const Rcpp::NumericVector& prob, const bool& verbose, const double& tol) {
  int x_int = static_cast<int>(x);
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
      if (s_k <= x_int) {
        Qk[0][s_k] = R::dbinom(s_k, size, prob[0], 0);
      }
    }
  }
  else if (k == m) {
    for (int s_km1 = 0; s_km1 <= size; s_km1++) {
      if (size - s_km1 <= x_int) {
        Qk[s_km1][0] = 1.0;
      }
    }
  }
  else {
    for (int s_km1 = 0; s_km1 <= size; s_km1++) {
      for (int s_k = s_km1; s_k <= size; s_k++) {
        if (s_k - s_km1 <= x_int) {
          Qk[s_km1][s_k] = R::dbinom(s_k - s_km1, size - s_km1, prob[k - 1]/prob_c[k - 1], 0);
        }
      }
    }
  }

  return Qk;
}

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
