#include <iostream>
#include <stdexcept>
#include <vector>

#include "JacobiSolver.hpp"
#include "SeidelSolver.hpp"

using namespace std;

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
    cout << "\nIterations: " << jacobi.getIterations() << endl;

    SeidelSolver seidel(A, b, eps);
    vector<double> x_seidel = seidel.solve();
    cout << "Seidel solution: ";
    for (double val : x_seidel) cout << val << " ";
    cout << "\nIterations: " << seidel.getIterations() << endl;

  } catch (const exception& e) {
    cerr << "Error: " << e.what() << endl;
  }
  return 0;
}