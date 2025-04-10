#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "newtons_method.cpp"
#include "simple_iteration.cpp"

double f(double x) { return tan(x) - 5 * x * x + 1; }

void summarize(const std::string& title, double x, int iterations) {
  std::cout << title << ":" << std::endl;
  std::cout << "x = " << std::fixed << std::setprecision(9) << x << std::endl;
  std::cout << "Iterations: " << iterations << std::endl;
  std::cout << "f(x) = " << f(x) << std::endl;
  std::cout << std::endl;
}

void print_errors(const std::string& method,
                  const std::vector<double>& errors) {
  std::cout << "Errors in " << method << ":" << std::endl;
  for (size_t i = 0; i < errors.size(); ++i) {
    std::cout << "Iter " << i + 1 << ": " << std::scientific << errors[i]
              << std::endl;
  }
  std::cout << std::endl;
}

int main() {
  double a = 0.4, b = 0.6;
  double initial_guess = (a + b) / 2;
  double epsilon = 1e-9;
  int max_iterations = 100;

  SimpleIterationsSolver simple(max_iterations);
  double x_simple = simple.solve(initial_guess, epsilon);
  int iters_simple = simple.get_iterations();
  summarize("Simple Iteration Method", x_simple, iters_simple);

  NewtonSolver newton(max_iterations);
  double x_newton = newton.solve(a, b, epsilon);
  int iters_newton = newton.get_iterations();
  summarize("Newton's Method", x_newton, iters_newton);

  auto errors_simple = simple.get_errors();
  print_errors("Simple Iteration", errors_simple);
  auto errors_newton = newton.get_errors();
  print_errors("Newton's Method", errors_newton);
}