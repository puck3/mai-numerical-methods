#pragma once
#include <stdexcept>
#include <vector>

class Matrix {
 private:
  std::vector<std::vector<double>> data;
  size_t rows;
  size_t cols;

 public:
  Matrix(size_t n = 0, size_t m = 0, double init = 0.0);
  Matrix(const std::vector<std::vector<double>>& d);

  size_t getRows() const;
  size_t getCols() const;
  double& operator()(size_t i, size_t j);
  const double& operator()(size_t i, size_t j) const;

  Matrix transpose() const;
  static Matrix identity(size_t n);
  void print(const std::string& title = "") const;

  Matrix operator*(const Matrix& other) const;
};