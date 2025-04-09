#include "QRDecomposer.hpp"

#include <cmath>

QRDecomposer::QRDecomposer(const Matrix& A) {
  R = A;
  Q = Matrix::identity(A.getRows());
  householder_transform(R);
}

double sign(double x) {
  double sign = x < 0 ? -1 : 1;
  return sign;
}

double subDiagColNorm(const Matrix& m, size_t j) {
  double norm = 0;
  for (size_t i = j; i < m.getRows(); ++i) {
    norm += m(i, j) * m(i, j);
  }
  return sqrt(norm);
}

std::vector<double> QRDecomposer::calculate_v(size_t j) const {
  const size_t n = R.getRows();
  std::vector<double> v(n, 0.0);

  v[j] = R(j, j) + sign(R(j, j)) * subDiagColNorm(R, j);
  for (size_t i = j + 1; i < n; ++i) {
    v[i] = R(i, j);
  }
  return v;
}

Matrix QRDecomposer::calculate_H(const std::vector<double>& v) {
  size_t n = v.size();
  Matrix H = Matrix::identity(n);
  double dot = 0;
  for (const auto& v_i : v) {
    dot += v_i * v_i;
  }
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      H(i, j) -= (2 * v[i] * v[j]) / dot;
    }
  }
  return H;
}

void QRDecomposer::householder_transform(Matrix& A) {
  const size_t n = A.getRows();

  for (size_t k = 0; k < n - 1; ++k) {
    auto v = calculate_v(k);
    Matrix H = calculate_H(v);
    R = H * R;
    Q = H * Q;
  }
}

Matrix QRDecomposer::getQ() const { return Q.transpose(); }
Matrix QRDecomposer::getR() const { return R; }