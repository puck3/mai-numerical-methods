#pragma once

#include <iostream>
#include <vector>

class Matrix {
 private:
  std::vector<std::vector<double>> data;
  size_t rows, cols;

 public:
  Matrix(size_t n = 0, size_t m = 0, double init = 0.0);

  Matrix(const std::vector<std::vector<double>>& d);

  size_t getRows() const;
  size_t getCols() const;

  double& operator()(size_t i, size_t j);
  const double& operator()(size_t i, size_t j) const;

  void swap_rows(size_t i, size_t j);

  std::vector<double>& operator[](size_t i);
  const std::vector<double>& operator[](size_t i) const;

  void print(const std::string& title = "") const;
};