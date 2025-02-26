#include <iostream>
#include <vector>

#include "LUDecomposer.hpp"
#include "Matrix.hpp"

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

    std::cout << "Determinant: " << lu.determinant() << std::endl;

    Matrix inv = lu.inverse();
    inv.print("Inverse matrix:");

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
}