#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

class SimpleIterationsSolver {
 private:
  std::vector<double> errors;
  int max_iterations;
  double g(double x) { return sqrt((tan(x) + 1) / 5.0); }

  // g'(x) = 1 / (2 * sqrt(5) * cos^2(x) * sqrt(tg(x) + 1))
  // g'(x)∈(0.2; 0.3), ∀x∈[0.4; 0.6]
  double q = 0.3;

 public:
  SimpleIterationsSolver(int max_iterations)
      : max_iterations(max_iterations), errors() {}

  double solve(double initial_guess, double epsilon) {
    double x = initial_guess;
    errors.clear();

    for (int i = 0; i < max_iterations; ++i) {
      double x_new = g(x);
      double error = (q * fabs(x_new - x)) / (1 - q);
      errors.push_back(error);
      if (error < epsilon) {
        return x_new;
      }
      x = x_new;
    }
    return x;
  }

  int get_iterations() const { return errors.size(); }

  const std::vector<double>& get_errors() const { return errors; }
};
