#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;

class Matrix {
 private:
  vector<vector<double>> data;
  size_t rows, cols;

 public:
  Matrix(size_t n = 0, size_t m = 0, double init = 0.0)
      : rows(n), cols(m), data(n, vector<double>(m, init)) {}

  Matrix(const vector<vector<double>>& d)
      : data(d), rows(d.size()), cols(d.empty() ? 0 : d[0].size()) {}

  size_t getRows() const { return rows; }
  size_t getCols() const { return cols; }

  double& operator()(size_t i, size_t j) { return data[i][j]; }
  const double& operator()(size_t i, size_t j) const { return data[i][j]; }

  void print(const string& title = "") const {
    if (!title.empty()) cout << title << endl;
    for (const auto& row : data) {
      for (double val : row) printf("%10.6f ", val);
      cout << endl;
    }
  }
};

class QRDecomposer {
 private:
  Matrix Q;
  Matrix R;

  void gram_schmidt(const Matrix& A) {
    size_t n = A.getRows();
    Q = Matrix(n, n);
    R = Matrix(n, n);

    for (size_t j = 0; j < n; ++j) {
      vector<double> v(n);
      for (size_t i = 0; i < n; ++i) v[i] = A(i, j);

      for (size_t k = 0; k < j; ++k) {
        R(k, j) = 0.0;
        for (size_t i = 0; i < n; ++i) R(k, j) += Q(i, k) * A(i, j);
        for (size_t i = 0; i < n; ++i) v[i] -= R(k, j) * Q(i, k);
      }

      double norm = 0.0;
      for (double x : v) norm += x * x;
      norm = sqrt(norm);

      R(j, j) = norm;
      for (size_t i = 0; i < n; ++i) Q(i, j) = v[i] / norm;
    }
  }

 public:
  QRDecomposer(const Matrix& A) { gram_schmidt(A); }

  Matrix getQ() const { return Q; }
  Matrix getR() const { return R; }
};

class QREigenSolver {
 private:
  Matrix A;
  double epsilon;
  int max_iterations;

  Matrix multiply(const Matrix& a, const Matrix& b) {
    size_t n = a.getRows();
    Matrix res(n, n);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        for (size_t k = 0; k < n; ++k) res(i, j) += a(i, k) * b(k, j);
    return res;
  }

  bool is_converged() {
    size_t n = A.getRows();
    for (size_t i = 1; i < n; ++i)
      for (size_t j = 0; j < i - 1; ++j)
        if (fabs(A(i, j)) > epsilon) return false;
    return true;
  }

 public:
  QREigenSolver(const Matrix& mat, double eps = 1e-6, int max_iter = 1000)
      : A(mat), epsilon(eps), max_iterations(max_iter) {}

  vector<double> solve() {
    size_t n = A.getRows();
    for (int iter = 0; iter < max_iterations; ++iter) {
      QRDecomposer qr(A);
      Matrix Q = qr.getQ();
      Matrix R = qr.getR();
      A = multiply(R, Q);
      if (is_converged()) break;
    }

    vector<double> eigenvalues(n);
    for (size_t i = 0; i < n; ++i) eigenvalues[i] = A(i, i);
    return eigenvalues;
  }
};

int main() {
  try {
    Matrix A({{6, 5, -6}, {4, -6, 9}, {-6, 6, 1}});

    QREigenSolver solver(A, 1e-6, 1000);
    vector<double> eigvals = solver.solve();

    cout << "Eigenvalues:\n";
    for (double val : eigvals) printf("%.6f ", val);
    cout << endl;

  } catch (const exception& e) {
    cerr << "Error: " << e.what() << endl;
  }
  return 0;
}