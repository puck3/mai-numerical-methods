#include <iostream>
#include <stdexcept>
#include <vector>

#include "JacobiSolver.hpp"
#include "SeidelSolver.hpp"

using namespace std;

void check_solution(const vector<vector<double>>& A,
                    const std::vector<double>& b,
                    const std::vector<double>& x) {
  std::cout << "A * x = (";
  for (size_t i = 0; i < A.size(); ++i) {
    double res = 0;
    for (size_t j = 0; j < A.size(); ++j) {
      res += (A[i][j] * x[j]);
    }
    std::cout << res << " ";
    if (res != b[i]) {
      // throw std::runtime_error("Wrong answer");
    }
    res = 0;
  }
  std::cout << ")^T" << std::endl;
}

int main() {
  try {
    size_t n;
    cin >> n;
    vector<vector<double>> A(n, vector<double>(n));
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        cin >> A[i][j];
      }
    }
    vector<double> b(n);
    for (size_t i = 0; i < n; ++i) {
      cin >> b[i];
    }
    double eps;
    cin >> eps;

    JacobiSolver jacobi(A, b, eps);
    vector<double> x_jacobi = jacobi.solve();
    cout << "Jacobi solution: ";
    for (double val : x_jacobi) cout << val << " ";
    cout << "\nIterations: " << jacobi.get_iterations() << endl;
    check_solution(A, b, x_jacobi);

    SeidelSolver seidel(A, b, eps);
    vector<double> x_seidel = seidel.solve();
    cout << "Seidel solution: ";
    for (double val : x_seidel) cout << val << " ";
    cout << "\nIterations: " << seidel.get_iterations() << endl;
    check_solution(A, b, x_seidel);

  } catch (const exception& e) {
    cerr << "Error: " << e.what() << endl;
  }
  return 0;
}