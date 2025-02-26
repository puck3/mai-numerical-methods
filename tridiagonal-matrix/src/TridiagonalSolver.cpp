#include "TridiagonalSolver.hpp"

#include <stdexcept>
#include <vector>

void TridiagonalSolver::validate() {
  if (mainDiagonal.size() != (size_t)n ||
      upperDiagonal.size() != (size_t)n - 1 ||
      lowerDiagonal.size() != (size_t)n - 1 || b.size() != (size_t)n) {
    throw std::invalid_argument("Invalid input sizes");
  }
}

TridiagonalSolver::TridiagonalSolver(const std::vector<double>& lower,
                                     const std::vector<double>& main,
                                     const std::vector<double>& upper,
                                     const std::vector<double>& rhs)
    : lowerDiagonal(lower),
      mainDiagonal(main),
      upperDiagonal(upper),
      b(rhs),
      n(main.size()) {
  validate();
}

void TridiagonalSolver::solve() {
  std::vector<double> alpha(n);
  std::vector<double> beta(n);

  alpha[0] = -upperDiagonal[0] / mainDiagonal[0];
  beta[0] = b[0] / mainDiagonal[0];

  for (int i = 1; i < n - 1; ++i) {
    double denom = mainDiagonal[i] + lowerDiagonal[i - 1] * alpha[i - 1];
    if (abs(denom) < 1e-10) {
      throw std::runtime_error("Matrix is singular or unstable");
    }
    alpha[i] = -upperDiagonal[i] / denom;
    beta[i] = (b[i] - lowerDiagonal[i - 1] * beta[i - 1]) / denom;
  }

  double denom = mainDiagonal[n - 1] + lowerDiagonal[n - 2] * alpha[n - 2];
  if (abs(denom) < 1e-10) {
    throw std::runtime_error("Matrix is singular or unstable");
  }
  beta[n - 1] = (b[n - 1] - lowerDiagonal[n - 2] * beta[n - 2]) / denom;

  x.resize(n);
  x[n - 1] = beta[n - 1];

  for (int i = n - 2; i >= 0; --i) {
    x[i] = alpha[i] * x[i + 1] + beta[i];
  }
}

const std::vector<double>& TridiagonalSolver::getSolution() const { return x; }