#pragma once
#include "Matrix.hpp"

class QRDecomposer {
 private:
  Matrix Q;
  Matrix R;

  void householder_transform(Matrix& A);
  std::vector<double> calculate_v(size_t j) const;
  static Matrix calculate_H(const std::vector<double>& v);

 public:
  QRDecomposer(const Matrix& A);
  Matrix getQ() const;
  Matrix getR() const;
};