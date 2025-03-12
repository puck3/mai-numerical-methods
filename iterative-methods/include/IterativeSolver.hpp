#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

using namespace std;

class IterativeSolver {
 protected:
  vector<vector<double>> A;
  vector<double> b;
  vector<double> x;
  double epsilon;
  int max_iterations;
  int iterations;

  virtual void prepare() {
    x.assign(A.size(), 0.0);
    for (size_t i = 0; i < A.size(); ++i) {
      double sum = 0.0;
      for (size_t j = 0; j < A[i].size(); ++j) {
        if (i != j) sum += fabs(A[i][j]);
      }
      if (fabs(A[i][i]) <= sum) {
        throw runtime_error("No diagonal dominance");
      }
      x[i] = b[i] / A[i][i];
    }
  }

  double norm(const vector<double>& v) {
    double max = 0.0;
    for (double val : v) {
      if (fabs(val) > max) max = fabs(val);
    }
    return max;
  }

 public:
  IterativeSolver(const vector<vector<double>>& matrix,
                  const vector<double>& rhs, double eps = 1e-6,
                  int max_iter = 1000)
      : A(matrix),
        b(rhs),
        epsilon(eps),
        max_iterations(max_iter),
        iterations(0) {}

  virtual vector<double> solve() = 0;

  int get_iterations() const { return iterations; }
};