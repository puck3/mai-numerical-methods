#include <iostream>
#include <vector>

#include "TridiagonalSolver.hpp"

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>,
           std::vector<double>>
input() {
  int n;
  std::cin >> n;
  std::vector<double> a(n - 1);
  std::vector<double> b(n);
  std::vector<double> c(n - 1);
  std::vector<double> d(n);

  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      std::cin >> a[i - 1];
    }
    std::cin >> b[i];
    if (i < (n - 1)) {
      std::cin >> c[i];
    }
    std::cin >> d[i];
  }
  return {a, b, c, d};
}

int main() {
  try {
    std::vector<double> a, b, c, d;
    std::tie(a, b, c, d) = input();
    TridiagonalSolver solver(a, b, c, d);
    solver.solve();

    std::vector<double> x = solver.getSolution();
    std::cout << "Solution: ";
    for (double xi : x) std::cout << xi << " ";
    std::cout << std::endl;

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
}