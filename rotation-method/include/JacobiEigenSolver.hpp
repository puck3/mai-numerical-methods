#include <cmath>
#include <stdexcept>
#include <vector>

using namespace std;

class JacobiEigenSolver {
 private:
  vector<vector<double>> A;
  vector<vector<double>> V;
  double eps;
  int max_iterations;
  int iterations;
  double error;

  void validate() {
    size_t n = A.size();
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < i; ++j) {
        if (abs(A[i][j] - A[j][i]) > 1e-9) {
          throw runtime_error("Matrix is not symmetric");
        }
      }
    }
  }

  void find_max_off_diagonal(int& p, int& q) {
    double max_val = 0.0;
    error = 0.0;
    size_t n = A.size();
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        if (abs(A[i][j]) > max_val) {
          max_val = abs(A[i][j]);
          p = i;
          q = j;
        }
        error += A[i][j] * A[i][j];
      }
    }
    error = sqrt(error);
  }

  void rotate(int i, int j) {
    double a_ii = A[i][i];
    double a_jj = A[j][j];
    double a_ij = A[i][j];

    double phi = 0.5 * atan2(2 * a_ij, a_jj - a_ii);
    double c = cos(phi);
    double s = sin(phi);

    A[i][i] = c * c * a_ii - 2 * s * c * a_ij + s * s * a_jj;
    A[j][j] = s * s * a_ii + 2 * s * c * a_ij + c * c * a_jj;
    A[i][j] = A[j][i] = 0.0;

    size_t n = A.size();
    for (int k = 0; k < n; ++k) {
      if (k != i && k != j) {
        double a_ki = A[k][i];
        double a_kj = A[k][j];
        A[k][i] = A[i][k] = c * a_ki - s * a_kj;
        A[k][j] = A[j][k] = s * a_ki + c * a_kj;
      }
    }

    for (int k = 0; k < n; ++k) {
      double v_ki = V[k][i];
      double v_kj = V[k][j];
      V[k][i] = c * v_ki - s * v_kj;
      V[k][j] = s * v_ki + c * v_kj;
    }
  }

 public:
  JacobiEigenSolver(const vector<vector<double>>& mat, double eps = 1e-6,
                    int max_iter = 1000)
      : A(mat), eps(eps), max_iterations(max_iter), iterations(0), error(0) {
    validate();
    size_t n = A.size();
    V.resize(n, vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
      V[i][i] = 1.0;
    }
  }

  void solve() {
    int i, j;
    for (iterations = 0; iterations < max_iterations; ++iterations) {
      find_max_off_diagonal(i, j);
      if (error < eps) break;
      rotate(i, j);
      iterations++;
    }
  }

  vector<double> eigenvalues() const {
    size_t n = A.size();
    vector<double> eig(n);
    for (size_t i = 0; i < n; ++i) {
      eig[i] = A[i][i];
    }
    return eig;
  }

  vector<vector<double>> get_eigenvectors() const { return V; }

  int get_iterations() const { return iterations; }
  double get_error() const { return error; }
};
