#include "LUDecomposer.hpp"

void LUDecomposer::forward_substitution(const std::vector<double>& pb,
                                        std::vector<double>& y) const {
  int n = LU.getRows();
  for (int i = 0; i < n; ++i) {
    y[i] = pb[i];
    for (int j = 0; j < i; ++j) y[i] -= LU(i, j) * y[j];
  }
}

void LUDecomposer::backward_substitution(const std::vector<double>& y,
                                         std::vector<double>& x) const {
  int n = LU.getRows();
  for (int i = n - 1; i >= 0; --i) {
    x[i] = y[i];
    for (int j = i + 1; j < n; ++j) x[i] -= LU(i, j) * x[j];
    x[i] /= LU(i, i);
  }
}

LUDecomposer::LUDecomposer(const Matrix& A)
    : LU(A), pivots(A.getRows()), num_swaps(0) {
  int n = A.getRows();
  for (int i = 0; i < n; ++i) pivots[i] = i;

  for (int k = 0; k < n; ++k) {
    int max_row = k;
    double max_val = abs(LU(k, k));
    for (int i = k + 1; i < n; ++i) {
      if (abs(LU(i, k)) > max_val) {
        max_val = abs(LU(i, k));
        max_row = i;
      }
    }

    if (max_row != k) {
      LU.swap_rows(k, max_row);
      std::swap(pivots[k], pivots[max_row]);
      num_swaps++;
    }

    if (LU(k, k) == 0) throw std::runtime_error("Matrix is singular");

    for (int i = k + 1; i < n; ++i) {
      LU(i, k) /= LU(k, k);
      for (int j = k + 1; j < n; ++j) {
        LU(i, j) -= LU(i, k) * LU(k, j);
      }
    }
  }
}

std::vector<double> LUDecomposer::solve(const std::vector<double>& b) const {
  int n = LU.getRows();
  std::vector<double> pb(n);
  for (int i = 0; i < n; ++i) pb[i] = b[pivots[i]];

  std::vector<double> y(n), x(n);
  forward_substitution(pb, y);
  backward_substitution(y, x);
  return x;
}

double LUDecomposer::determinant() const {
  double det = 1.0;
  for (int i = 0; i < LU.getRows(); ++i) det *= LU(i, i);
  return det * (num_swaps % 2 ? -1 : 1);
}

Matrix LUDecomposer::inverse() const {
  int n = LU.getRows();
  Matrix inv(n, n);

  for (int col = 0; col < n; ++col) {
    std::vector<double> e(n, 0.0);
    e[col] = 1.0;

    std::vector<double> pb(n);
    for (int i = 0; i < n; ++i) pb[i] = e[pivots[i]];

    std::vector<double> y(n), x(n);
    forward_substitution(pb, y);
    backward_substitution(y, x);

    for (int row = 0; row < n; ++row) inv(row, col) = x[row];
  }
  return inv;
}
