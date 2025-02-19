#pragma once

#include <iostream>
#include <vector>

class SquareMatrix {
 private:
  size_t n;
  std::vector<std::vector<double>> data;

 public:
  SquareMatrix() = default;

  SquareMatrix(size_t n);

  virtual ~SquareMatrix() = default;

  size_t size() const;

  double& operator()(size_t i, size_t j);

  const double& operator()(size_t i, size_t j) const;

  friend std::ostream& operator<<(std::ostream& os, const SquareMatrix& m);
  friend std::istream& operator>>(std::istream& is, SquareMatrix& m);
};