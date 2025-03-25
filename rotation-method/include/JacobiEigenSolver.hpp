#include <cmath>
#include <stdexcept>
#include <vector>

using namespace std;

class JacobiEigenSolver {
 private:
  vector<vector<double>> matrix;
  vector<vector<double>> eigenvectors;
  double epsilon;
  int max_iterations;
  int iterations;
  double error;

  void validate() {
    size_t n = matrix.size();
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < i; ++j) {
        if (abs(matrix[i][j] - matrix[j][i]) > 1e-9) {
          throw runtime_error("Matrix is not symmetric");
        }
      }
    }
  }

  void find_max_off_diagonal(int& p, int& q) {
    double max_val = 0.0;
    size_t n = matrix.size();
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        if (abs(matrix[i][j]) > max_val) {
          max_val = abs(matrix[i][j]);
          p = i;
          q = j;
        }
      }
    }
    error = max_val;
  }

  void rotate(int p, int q) {
    double app = matrix[p][p];
    double aqq = matrix[q][q];
    double apq = matrix[p][q];

    double theta = 0.5 * atan2(2 * apq, aqq - app);
    double c = cos(theta);
    double s = sin(theta);

    matrix[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq;
    matrix[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq;
    matrix[p][q] = matrix[q][p] = 0.0;

    size_t n = matrix.size();
    for (int i = 0; i < n; ++i) {
      if (i != p && i != q) {
        double aip = matrix[i][p];
        double aiq = matrix[i][q];
        matrix[i][p] = matrix[p][i] = c * aip - s * aiq;
        matrix[i][q] = matrix[q][i] = s * aip + c * aiq;
      }
    }

    for (int i = 0; i < n; ++i) {
      double vip = eigenvectors[i][p];
      double viq = eigenvectors[i][q];
      eigenvectors[i][p] = c * vip - s * viq;
      eigenvectors[i][q] = s * vip + c * viq;
    }
  }

 public:
  JacobiEigenSolver(const vector<vector<double>>& mat, double eps = 1e-6,
                    int max_iter = 1000)
      : matrix(mat),
        epsilon(eps),
        max_iterations(max_iter),
        iterations(0),
        error(0) {
    validate();
    size_t n = matrix.size();
    eigenvectors.resize(n, vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
      eigenvectors[i][i] = 1.0;
    }
  }

  void solve() {
    int p, q;
    iterations = 0;
    do {
      find_max_off_diagonal(p, q);
      if (error < epsilon) break;
      rotate(p, q);
      iterations++;
    } while (iterations < max_iterations);
  }

  vector<double> eigenvalues() const {
    size_t n = matrix.size();
    vector<double> eig(n);
    for (size_t i = 0; i < n; ++i) {
      eig[i] = matrix[i][i];
    }
    return eig;
  }

  vector<vector<double>> get_eigenvectors() const { return eigenvectors; }

  int get_iterations() const { return iterations; }
  double get_error() const { return error; }
};
