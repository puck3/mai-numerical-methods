#include "TridiagonalSolver.hpp"

#include <stdexcept>
#include <vector>

void TridiagonalSolver::validate() {
  if (b.size() != (size_t)n || c.size() != (size_t)n - 1 ||
      a.size() != (size_t)n - 1 || d.size() != (size_t)n) {
    throw std::invalid_argument("Invalid input sizes");
  }
}

TridiagonalSolver::TridiagonalSolver(const std::vector<double>& a,
                                     const std::vector<double>& b,
                                     const std::vector<double>& c,
                                     const std::vector<double>& d)
    : a(a), b(b), c(c), d(d), n(b.size()) {
  validate();
}

void TridiagonalSolver::solve() {
  std::vector<double> P(n);
  std::vector<double> Q(n);

  P[0] = -c[0] / b[0];
  Q[0] = d[0] / b[0];

  for (int i = 1; i < n - 1; ++i) {
    double denom = b[i] + a[i - 1] * P[i - 1];
    if (abs(denom) < 1e-10) {
      throw std::runtime_error("Matrix is singular or unstable");
    }
    P[i] = -c[i] / denom;
    Q[i] = (d[i] - a[i - 1] * Q[i - 1]) / denom;
  }

  double denom = b[n - 1] + a[n - 2] * P[n - 2];
  if (abs(denom) < 1e-10) {
    throw std::runtime_error("Matrix is singular or unstable");
  }
  Q[n - 1] = (d[n - 1] - a[n - 2] * Q[n - 2]) / denom;

  x.resize(n);
  x[n - 1] = Q[n - 1];

  for (int i = n - 2; i >= 0; --i) {
    x[i] = P[i] * x[i + 1] + Q[i];
  }
}

const std::vector<double>& TridiagonalSolver::getSolution() const { return x; }