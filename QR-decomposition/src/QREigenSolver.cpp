#include "QREigenSolver.hpp"

#include <cmath>

#include "QRDecomposer.hpp"

QREigenSolver::QREigenSolver(const Matrix& mat, double eps, int max_iter)
    : A(mat), epsilon(eps), max_iterations(max_iter) {}

bool QREigenSolver::is_converged() const {
  const size_t n = A.getRows();
  double epsilon_k = 0.0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j + 1 < i; ++j) {
      epsilon_k += A(i, j) * A(i, j);
    }
  }
  epsilon_k = sqrt(epsilon_k);
  if (epsilon_k <= epsilon) {
    return true;
  } else {
    return false;
  }
}

std::vector<double> QREigenSolver::solve() {
  const size_t n = A.getRows();

  for (int iter = 0; iter < max_iterations; ++iter) {
    QRDecomposer qr(A);
    Matrix Q = qr.getQ();
    Matrix R = qr.getR();
    A = R * Q;
    if (is_converged()) {
      break;
    }
  }

  std::vector<double> eigenvalues(n);
  for (size_t i = 0; i < n; ++i) {
    eigenvalues[i] = A(i, i);
  }
  return eigenvalues;
}
