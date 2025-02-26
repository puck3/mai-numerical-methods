#pragma once

#include <cmath>
#include <vector>

#include "Matrix.hpp"

class LUDecomposer {
 private:
  Matrix LU;
  std::vector<int> pivots;
  int num_swaps;

  void forward_substitution(const std::vector<double>& pb,
                            std::vector<double>& y) const;

  void backward_substitution(const std::vector<double>& y,
                             std::vector<double>& x) const;

 public:
  LUDecomposer(const Matrix& A);

  std::vector<double> solve(const std::vector<double>& b) const;

  double determinant() const;

  Matrix inverse() const;
};
