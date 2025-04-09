#pragma once
#include <vector>

#include "Matrix.hpp"

class QREigenSolver {
 private:
  Matrix A;
  double epsilon;
  int max_iterations;

  bool is_converged() const;

 public:
  QREigenSolver(const Matrix& mat, double eps = 1e-6, int max_iter = 1000);

  std::vector<double> solve();
};