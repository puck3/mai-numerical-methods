#pragma once

#include <iostream>
#include <vector>

class Matrix {
 private:
  size_t rows;
  size_t cols;
  std::vector<std::vector<double>> data;

 public:
  Matrix() = default;

  Matrix(size_t n, size_t m);

  virtual ~Matrix() = default;

  double& operator()(size_t i, size_t j);

  const double& operator()(size_t i, size_t j) const;

  friend std::ostream& operator<<(std::ostream& os, const Matrix& m);
  friend std::istream& operator>>(std::istream& is, Matrix& m);
};