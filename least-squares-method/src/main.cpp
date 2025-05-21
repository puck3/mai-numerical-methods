#include <cmath>
#include <iostream>
#include <vector>

#include "LUDecomposer.hpp"
#include "Matrix.hpp"

class LeastSquares {
 private:
  std::vector<double> a;
  size_t n;

 public:
  LeastSquares(const std::vector<double>& x, const std::vector<double>& y,
               size_t n)
      : n(n) {
    std::vector<std::vector<double>> MatrixData(n + 1,
                                                std::vector<double>(n + 1, 0));
    std::vector<double> b(n + 1, 0);
    for (size_t i = 0; i <= n; ++i) {
      for (size_t j = 0; j <= n; ++j) {
        for (size_t k = 0; k < x.size(); ++k) {
          MatrixData[i][j] += pow(x[k], i + j);
        }
      }
      for (size_t k = 0; k < x.size(); ++k) {
        b[i] += y[k] * pow(x[k], i);
      }
    }
    Matrix A(MatrixData);
    LUDecomposer lu(A);
    a = lu.solve(b);
  }

  void print() const {
    for (size_t i = 0; i <= n; ++i) {
      std::cout << a[i] << "x^" << i;
      if (i != n) {
        std::cout << " + ";
      }
    }
    std::cout << std::endl;
  }

  double operator()(double x) const {
    double res = 0;
    for (size_t i = 0; i <= n; ++i) {
      res += a[i] * pow(x, i);
    }
    return res;
  }
};

int main() {
  try {
    size_t n;
    std::cin >> n;
    std::vector<double> x(n);
    std::vector<double> y(n);
    for (size_t i = 0; i < n; ++i) {
      std::cin >> x[i];
    }
    for (size_t i = 0; i < n; ++i) {
      std::cin >> y[i];
    }

    size_t q;
    std::cin >> q;
    for (size_t _ = 0; _ < q; ++_) {
      size_t m;
      std::cin >> m;
      LeastSquares f(x, y, m);
      f.print();
      double square_error_sum = 0;
      for (size_t i = 0; i < n; ++i) {
        square_error_sum += pow(f(x[i]) - y[i], 2);
      }
      std::cout << "Error: " << square_error_sum << std::endl;
    }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
}