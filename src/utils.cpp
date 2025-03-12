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

std::vector<std::vector<double>> addMatrix(const std::vector<std::vector<double>>& A,
  const std::vector<std::vector<double>>& B) {
  int rowsA = A.size();
  int colsA = A[0].size();
  int rowsB = B.size();
  int colsB = B[0].size();

  // ensure matrix addition is valid
  if (rowsA != rowsB) {
    throw std::invalid_argument("number of rows of A must match number of rows of B.");
  }
  if (colsA != colsB) {
    throw std::invalid_argument("number of columns of A must match number of columns of B.");
  }

  // initialize C matrix with zeroes
  std::vector<std::vector<double>> C(rowsA, std::vector<double>(colsA, 0.0));

  // perform matrix multiplication
  for (int i = 0; i < rowsA; i++) {
    for (int j = 0; j < colsA; j++) {
      C[i][j] = A[i][j] + B[i][j];
    }
  }

  return C;
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

// [[Rcpp::export]]
Rcpp::NumericMatrix vector2D_2_NM(std::vector<std::vector<double>> mat) {
  int rows = mat.size();
  if (rows == 0) {
    Rcpp::stop("input matrix has no rows.");
  }
  
  int cols = mat[0].size();
  for (const auto& row : mat) {
    if (row.size() != cols) {
      Rcpp::stop("all rows must have the same number of columns.");
    }
  }

  Rcpp::NumericMatrix result(rows, cols);  // create an empty Rcpp::NumericMatrix

  // fill NumericMatrix by converting from row-major to column-major format
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      result(i, j) = mat[i][j];  // direct assignment
    }
  }

  return result;
}

// [[Rcpp::export]]
arma::mat vec_2_armaMat(std::vector<std::vector<double>> x) {
  // Determine matrix dimensions
  size_t rows = x.size();
  size_t cols = rows > 0 ? x[0].size() : 0;

  // Create an Armadillo matrix
  arma::mat mat(rows, cols);

  // Fill the matrix
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      mat(i, j) = x[i][j];
    }
  }

  return mat;
}

// [[Rcpp::export]]
std::vector<std::vector<double>> armaMat_2_vec(const arma::mat& mat) {
  size_t rows = mat.n_rows;
  size_t cols = mat.n_cols;

  // Create a nested vector with the same dimensions
  std::vector<std::vector<double>> vec(rows, std::vector<double>(cols));

  // Copy data from arma::mat to std::vector<std::vector<double>>
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      vec[i][j] = mat(i, j);
    }
  }

  return vec;
}
