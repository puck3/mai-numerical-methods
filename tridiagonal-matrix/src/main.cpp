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

void check_solution(const std::vector<double>& a, const std::vector<double>& b,
                    const std::vector<double>& c, const std::vector<double>& d,
                    const std::vector<double>& x) {
  std::cout << "A * x = (";
  double di = 0;
  for (int i = 0; i < x.size(); ++i) {
    if (i > 0) {
      di += x[i - 1] * a[i - 1];
    }
    di += x[i] * b[i];
    if (i < x.size() - 1) {
      di += x[i + 1] * c[i];
    }
    if (di != d[i]) {
      throw std::runtime_error("Wrong answer");
    }
    std::cout << di << " ";
    di = 0;
  }
  std::cout << ")^T" << std::endl;
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

    check_solution(a, b, c, d, x);

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
}