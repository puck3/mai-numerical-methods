#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "NonlinearSolver.hpp"

double f1(double x1, double x2) { return x1 * x1 - 2 * log10(x2) - 1; }

double f2(double x1, double x2) { return x1 * x1 - 2 * x1 * x2 + 2; }

int main() {
  try {
    NonlinearSolver solver(1e-9, 100);

    double x1_init = 1.0;
    double x2_init = 1.5;

    Result si_result = solver.simpleIteration(x1_init, x2_init);
    Result nm_result = solver.newtonMethod(x1_init, x2_init);

    std::cout << "Simple Iteration Method:" << std::endl;
    std::cout << "Iterations: " << si_result.iterations << std::endl;
    std::cout << "Solution: x1 = " << si_result.x1 << ", x2 = " << si_result.x2
              << std::endl;
    std::cout << "f_1(x1, x2) = " << f1(si_result.x1, si_result.x2)
              << std::endl;
    std::cout << "f_2(x1, x2) = " << f2(si_result.x1, si_result.x2) << std::endl
              << std::endl;

    std::cout << "Newton's Method:" << std::endl;
    std::cout << "Iterations: " << nm_result.iterations << std::endl;
    std::cout << "Solution: x1 = " << nm_result.x1 << ", x2 = " << nm_result.x2
              << std::endl;
    std::cout << "f_1(x1, x2) = " << f1(nm_result.x1, nm_result.x2)
              << std::endl;
    std::cout << "f_2(x1, x2) = " << f2(nm_result.x1, nm_result.x2) << std::endl
              << std::endl;

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}