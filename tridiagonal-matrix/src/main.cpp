#include <iostream>
#include <vector>

#include "TridiagonalSolver.hpp"

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>,
           std::vector<double>>
input() {
  int n;
  std::cin >> n;
  std::vector<double> lowerDiagonal(n - 1);
  std::vector<double> mainDiagonal(n);
  std::vector<double> upperDiagonal(n - 1);
  std::vector<double> b(n);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double elem;
      std::cin >> elem;
      if (i == j) {
        mainDiagonal[i] = elem;
      } else if (i == j + 1) {
        lowerDiagonal[j] = elem;
      } else if (i + 1 == j) {
        upperDiagonal[i] = elem;
      } else if (elem != 0) {
        throw std::runtime_error("Matrix is not tridiagonal");
      }
    }
    std::cin >> b[i];
  }
  return {lowerDiagonal, mainDiagonal, upperDiagonal, b};
}

int main() {
  try {
    std::vector<double> lowerDiagonal, mainDiagonal, upperDiagonal, b;
    std::tie(lowerDiagonal, mainDiagonal, upperDiagonal, b) = input();
    TridiagonalSolver solver(lowerDiagonal, mainDiagonal, upperDiagonal, b);
    solver.solve();

    std::vector<double> x = solver.getSolution();
    std::cout << "Solution: ";
    for (double xi : x) std::cout << xi << " ";
    std::cout << std::endl;

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
}