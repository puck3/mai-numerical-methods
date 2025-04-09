#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "JacobiEigenSolver.hpp"

using namespace std;

vector<double> matrix_vector_mult(const vector<vector<double>>& A,
                                  const vector<double>& v) {
  size_t n = A.size();
  vector<double> result(n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      result[i] += A[i][j] * v[j];
    }
  }
  return result;
}

vector<double> scalar_vector_mult(double a, const vector<double>& v) {
  size_t n = v.size();
  vector<double> result(n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    result[i] = a * v[i];
  }
  return result;
}

vector<vector<double>> transpose(const vector<vector<double>>& mat) {
  size_t n = mat.size();
  vector<vector<double>> transposed(n, vector<double>(n));
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      transposed[j][i] = mat[i][j];
    }
  }
  return transposed;
}

void test_eigen_vectors_and_values(vector<vector<double>> A,
                                   vector<vector<double>> vecs,
                                   vector<double> vals) {
  vecs = transpose(vecs);
  for (int i = 0; i < vecs.size(); ++i) {
    vector<double> res1 = matrix_vector_mult(A, vecs[i]);
    cout << "A * v_" << i << " = (";
    for (int j = 0; j < res1.size(); ++j) {
      cout << res1[j];
      if (j < res1.size() - 1) {
        cout << " ";
      } else {
        cout << ")^T" << endl;
      }
    }

    vector<double> res2 = scalar_vector_mult(vals[i], vecs[i]);
    cout << "λ_" << i << " * v_" << i << " = (";
    for (int j = 0; j < res2.size(); ++j) {
      cout << res2[j];
      if (j < res1.size() - 1) {
        cout << " ";
      } else {
        cout << ")^T" << endl;
      }
    }
    cout << endl;
  }
}

int main() {
  try {
    int n;
    cin >> n;
    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        cin >> A[i][j];
      }
    }

    double eps;
    cin >> eps;

    JacobiEigenSolver solver(A, eps);
    solver.solve();

    vector<double> eigvals = solver.eigenvalues();
    vector<vector<double>> eigvecs = solver.get_eigenvectors();

    cout << "Eigenvalues:\n";
    for (double val : eigvals) cout << val << " ";
    cout << "\n\nEigenvectors:\n";
    for (auto& vec : eigvecs) {
      for (double v : vec) cout << v << "\t";
      cout << endl;
    }

    cout << "\nIterations: " << solver.get_iterations();
    cout << "\nFinal error: " << solver.get_error() << endl;

    test_eigen_vectors_and_values(A, eigvecs, eigvals);

  } catch (const exception& e) {
    cerr << "Error: " << e.what() << endl;
  }
  return 0;
}