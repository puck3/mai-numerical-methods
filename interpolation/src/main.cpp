#include <cmath>
#include <iostream>
#include <vector>

#include "LagrangePolynomial.hpp"
#include "NewtonPolynomial.hpp"

double y(double x) { return acos(x) + x; }

int main() {
  double n;
  std::cin >> n;
  std::vector<double> X_a(n);
  for (size_t i = 0; i < n; i++) {
    std::cin >> X_a[i];
  }
  std::vector<double> X_b(n);
  for (size_t i = 0; i < n; i++) {
    std::cin >> X_b[i];
  }
  double x;
  std::cin >> x;

  LagrangePolynomial lagrange(y, X_a);
  NewtonPolynomial newton(y, X_b);
  double res;

  std::cout << "a) X = ";
  for (double x_i : X_a) {
    std::cout << x_i << " ";
  }
  std::cout << std::endl;

  res = lagrange.interpolate(x);
  std::cout << "Lagrange Polynomial:" << std::endl;
  std::cout << "P(x*) = " << res << std::endl;
  std::cout << "Error: " << fabs(y(x) - res) << std::endl;

  res = newton.interpolate(x);
  std::cout << "Newton Polynomial:" << std::endl;
  std::cout << "P(x*) = " << res << std::endl;
  std::cout << "Error: " << fabs(y(x) - res) << std::endl;
  std::cout << std::endl;

  LagrangePolynomial lagrange2(y, X_b);
  NewtonPolynomial newton2(y, X_a);
  std::cout << "b) X = ";
  for (double x_i : X_b) {
    std::cout << x_i << " ";
  }
  std::cout << std::endl;

  res = lagrange2.interpolate(x);
  std::cout << "Lagrange Polynomial:" << std::endl;
  std::cout << "P(x*) = " << res << std::endl;
  std::cout << "Error: " << fabs(y(x) - res) << std::endl;

  res = newton2.interpolate(x);
  std::cout << "Newton Polynomial:" << std::endl;
  std::cout << "P(x*) = " << res << std::endl;
  std::cout << "Error: " << fabs(y(x) - res) << std::endl;
  std::cout << std::endl;
  return 0;
}