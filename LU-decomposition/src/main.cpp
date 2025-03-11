#include <iostream>
#include <vector>

#include "LUDecomposer.hpp"
#include "Matrix.hpp"

const int eps = 1e-9;

void check_solution(const Matrix& A, const std::vector<double>& b,
                    const std::vector<double>& x) {
  std::cout << "A * x = (";
  for (size_t i = 0; i < A.getRows(); ++i) {
    double res = 0;
    for (size_t j = 0; j < A.getRows(); ++j) {
      res += (A[i][j] * x[j]);
    }
    std::cout << res << " ";
    if (res != b[i]) {
      // throw std::runtime_error("Invalid answer");
    }
    res = 0;
  }
  std::cout << ")^T" << std::endl;
}

void check_inverse(const Matrix& A, const Matrix& InversedA) {
  double res = 0;
  std::cout << "A * A^(-1) = " << std::endl;
  for (size_t i = 0; i < A.getRows(); ++i) {
    for (size_t j = 0; j < A.getCols(); ++j) {
      for (size_t k = 0; k < A.getRows(); ++k) {
        res += A[i][k] * InversedA[k][j];
      }
      if ((i == j && abs(res - 1) > eps) || (i != j && abs(res) > eps)) {
        throw std::runtime_error("Invalid inversed matrix");
      }
      std::cout << res << " ";
      res = 0;
    }
    std::cout << std::endl;
  }
}

int main() {
  try {
    size_t n;
    std::cin >> n;
    std::vector<std::vector<double>> MatrixData(n, std::vector<double>(n));
    std::vector<double> b(n);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        std::cin >> MatrixData[i][j];
      }
      std::cin >> b[i];
    }
    Matrix A(MatrixData);
    LUDecomposer lu(A);

    std::vector<double> x = lu.solve(b);
    std::cout << "Solution x: ";
    for (double xi : x) std::cout << xi << " ";
    std::cout << std::endl;
    check_solution(A, b, x);

    std::cout << "Determinant: " << lu.determinant() << std::endl;

    Matrix inv = lu.inverse();
    inv.print("Inverse matrix:");
    check_inverse(A, inv);

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
}