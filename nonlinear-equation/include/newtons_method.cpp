#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

class NewtonSolver {
 private:
  const double max_iterations;
  std::vector<double> errors;

  double f(double x) { return tan(x) - 5 * x * x + 1; }

  double g(double x) { return (1.0 / (cos(x) * cos(x))) - 10 * x; }

  double h(double x) { return (2 * sin(x)) / (cos(x) * cos(x) * cos(x)) - 10; }

 public:
  NewtonSolver(double max_iterations)
      : max_iterations(max_iterations), errors() {}

  double solve(double a, double b, double eps) {
    double x;
    if (f(a) * f(b) > 0) {
      throw std::runtime_error("");
    }
    if (f(a) * h(a) > 0) {
      x = a;
    } else if (f(b) * h(b) > 0) {
      x = b;
    } else {
      throw std::runtime_error("Invalid value");
    }
    errors.clear();

    for (int i = 0; i < max_iterations; ++i) {
      double fx = f(x);
      double dfx = g(x);

      if (fabs(dfx) < 1e-12) {
        throw std::runtime_error("Derivative is too small");
      }

      double x_new = x - fx / dfx;
      double error = fabs(x_new - x);
      errors.push_back(error);

      if (error < eps) {
        return x_new;
      }
      x = x_new;
    }
    return x;
  }

  const std::vector<double>& get_errors() const { return errors; }

  int get_iterations() const { return errors.size(); }
};
