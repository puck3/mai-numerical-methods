#include <iostream>
#include <stdexcept>
#include <vector>

#include "Matrix.hpp"
#include "QREigenSolver.hpp"

Matrix readMatrix() {
  int n;
  std::cin >> n;
  std::vector<std::vector<double>> m(n, std::vector<double>(n));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      std::cin >> m[i][j];
    }
  }
  return Matrix(m);
}

double readEps() {
  double eps;
  std::cin >> eps;
  return eps;
}

void test_eigenvalues(const std::vector<double>& vals, const Matrix& A) {
  size_t n = A.getRows();
  if (n == 3) {
    double res = 1;
    for (size_t i = 0; i < n; ++i) {
      std::cout << "λ_" << i << (i == n - 1 ? " = " : " * ");
      res *= vals[i];
    }
    std::cout << res << std::endl;
    double det = A(0, 0) * A(1, 1) * A(2, 2);
    det += A(0, 1) * A(1, 2) * A(2, 0);
    det += A(0, 2) * A(1, 0) * A(2, 1);
    det -= A(0, 2) * A(1, 1) * A(2, 0);
    det -= A(0, 1) * A(1, 0) * A(2, 2);
    det -= A(0, 0) * A(1, 2) * A(2, 1);
    std::cout << "det(A) = " << det << std::endl;
  }

  double sum = 0, track = 0;
  for (size_t i = 0; i < n; ++i) {
    track += A(i, i);
    sum += vals[i];
    std::cout << "λ_" << i << (i == n - 1 ? " = " : " + ");
  }
  std::cout << sum << std::endl;
  std::cout << "tr(A) = " << track << std::endl;
}

int main() {
  try {
    Matrix A = readMatrix();
    double eps = readEps();

    QREigenSolver solver(A, eps, 1000);
    auto eigenvalues = solver.solve();

    std::cout << "Eigenvalues:" << std::endl;
    for (double val : eigenvalues) {
      std::cout << val << " ";
    }
    std::cout << std::endl;

    test_eigenvalues(eigenvalues, A);
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  return 0;
}