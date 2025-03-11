// -*- mode: C++; c-indent-level: 2; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "XOMultinom.h"

// Note: RcppExport is an alias for extern "C"

// [[Rcpp::export]]
void twoloops(const int& d, const int& n) {
  std::vector<int> s(d);
  std::vector<int> np(d + 1);
  int k = 1;

  np[0] = 1;
  for (int u = 1; u <= d; ++u) {
    np[u] = np[u - 1] * n;
    k *= n;
  }

  int tmp = 0;
  for (int u = 0; u < k; ++u) {
    for (int v = 0; v < d - 1; ++v) {
      tmp = 1 + (u % np[v + 1]) / np[v];
      s[v] = tmp;
      // std::cout << s[v] << " ";
    }
    // std::cout << std::endl;

    R_CheckUserInterrupt();
  }
}

// [[Rcpp::export]]
void twoloops_matrix(const int& d, const int& n) {
  std::vector<int> s(d);
  std::vector<int> np(d + 1);

  for (int u = 0; u <= d; ++u) {
    np[u] = std::pow(n, u);
  }

  int k = np[d];

  // Prepare the transformation matrix manually
  std::vector<std::vector<int>> M(d - 1, std::vector<int>(d, 0));
  for (int v = 0; v < d - 1; ++v) {
    M[v][v] = 1;
    M[v][v + 1] = -1;
  }

  for (int u = 0; u < k; ++u) {
    std::vector<int> U(d);
    for (int i = 0; i < d; i++) {
      U[i] = (u % np[i + 1]) / np[i];
    }

    // matrix-vector multiplication
    for (int v = 0; v < d - 1; ++v) {
      int result = 0;
      for (int j = 0; j < d; j++) {
        result += M[v][j] * U[j];
      }
      s[v] = 1 + result;
      // std::cout << s[v] << " ";
    }
    // std::cout << std::endl;

    R_CheckUserInterrupt();
  }
}

void action(std::vector<double> indices) {
  for (int i = 0; i < indices.size(); ++i) {
    std::cout << indices[i] << " ";
  }
  std::cout << std::endl;
}

void create_loops(int current_level, std::vector<double> indices, double xa, int levels,
  void (*act)(std::vector<double>)) {
  if (current_level > levels) {
    // base case: all loops are complete, perform the action
    act(indices);
  } else {
    // recursive case: create the current loop and recurse
    for (int i = 0; i <= xa; i++) {
      indices[current_level - 1] = i;
      create_loops(current_level + 1, indices, xa, levels, act);

      R_CheckUserInterrupt();
    }
  }
}

// [[Rcpp::export]]
void dynamic_nested_loops(int levels, double x) {
  int one = 1;
  std::vector<double> indices(levels, -999);
  create_loops(one, indices, x, levels, action);
}
