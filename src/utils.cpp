// -*- mode: C++; c-indent-level: 2; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "XOMultinom.h"

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
Rcpp::DataFrame aggregate_sum_rcpp(Rcpp::DataFrame x, Rcpp::List by) {
  Rcpp::Function f("aggregate_sum_df");

  return f(Rcpp::Named("x") = x, Rcpp::Named("by") = by);
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

void print_umap(std::unordered_map<int, double> myMap) {
  // Iterating through the hash table
  for (const auto& pair : myMap) {
    std::cout << pair.first << " -> " << pair.second << "\n";
  }
}

std::vector<std::vector<double>> multiplyMatrix(const std::vector<std::vector<double>>& A,
  const std::vector<std::vector<double>>& B) {
  int rowsA = A.size();
  int colsA = A[0].size();
  int rowsB = B.size();
  int colsB = B[0].size();

  // ensure matrix multiplication is valid
  if (colsA != rowsB) {
    throw std::invalid_argument("number of columns of A must match number of rows of B.");
  }

  // initialize C matrix with zeroes
  std::vector<std::vector<double>> C(rowsA, std::vector<double>(colsB, 0.0));

  // perform matrix multiplication
  for (int i = 0; i < rowsA; i++) {
    for (int j = 0; j < colsB; j++) {
      for (int k = 0; k < colsA; k++) {
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }

  return C;
}

std::vector<std::vector<double>> multiplyUpperTriangular(const std::vector<std::vector<double>>& A,
  const std::vector<std::vector<double>>& B) {
  int rowsA = A.size();
  int colsA = A[0].size();
  int rowsB = B.size();
  int colsB = B[0].size();

  // ensure matrix multiplication is valid
  if (colsA != rowsB) {
    throw std::invalid_argument("number of columns of A must match number of rows of B.");
  }
  if (colsA != rowsA) {
    throw std::invalid_argument("matrix A must be square.");
  }
  if (colsB != rowsB) {
    throw std::invalid_argument("matrix B must be square.");
  }

  // initialize C matrix with zeroes
  std::vector<std::vector<double>> C(rowsA, std::vector<double>(rowsA, 0.0));

  for (int i = 0; i < rowsA; ++i) {
    for (int j = i; j < rowsA; ++j) {  // j starts from i (upper triangular property)
      C[i][j] = 0;
      for (int k = i; k <= j; ++k) {  // k starts from i
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }

  return C;
}

// Function to print a matrix
void printMatrix(const std::vector<std::vector<double>>& matrix) {
  for (const auto& row : matrix) {
    for (double val : row) {
      std::cout << val << " ";
    }
    std::cout << std::endl;
  }
}
