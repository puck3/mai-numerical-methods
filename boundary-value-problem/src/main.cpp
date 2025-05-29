#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

double exact_y(double x) { return 2 * x + 1 + exp(2 * x); }

void derivs(double x, const vector<double>& y, vector<double>& dydx) {
  dydx.resize(2);
  if (fabs(x) < 1e-9) {
    dydx[0] = y[1];
    dydx[1] = (2 * x + 1) * y[1] - 2 * y[0];
  } else {
    dydx[0] = y[1];
    dydx[1] = ((2 * x + 1) * y[1] - 2 * y[0]) / x;
  }
}

void rk4_system(double x0, vector<double> y0, double xn, double h,
                vector<double>& x_vec, vector<vector<double>>& y_vec) {
  x_vec.clear();
  y_vec.clear();

  double x = x0;
  vector<double> y = y0;

  x_vec.push_back(x);
  y_vec.push_back(y);

  int steps = static_cast<int>((xn - x0) / h + 1e-9);

  for (int i = 0; i < steps; ++i) {
    vector<double> dydx;
    derivs(x, y, dydx);

    vector<double> k1(2);
    k1[0] = h * dydx[0];
    k1[1] = h * dydx[1];

    vector<double> y_temp(2);
    y_temp[0] = y[0] + 0.5 * k1[0];
    y_temp[1] = y[1] + 0.5 * k1[1];
    derivs(x + 0.5 * h, y_temp, dydx);

    vector<double> k2(2);
    k2[0] = h * dydx[0];
    k2[1] = h * dydx[1];

    y_temp[0] = y[0] + 0.5 * k2[0];
    y_temp[1] = y[1] + 0.5 * k2[1];
    derivs(x + 0.5 * h, y_temp, dydx);

    vector<double> k3(2);
    k3[0] = h * dydx[0];
    k3[1] = h * dydx[1];

    y_temp[0] = y[0] + k3[0];
    y_temp[1] = y[1] + k3[1];
    derivs(x + h, y_temp, dydx);

    vector<double> k4(2);
    k4[0] = h * dydx[0];
    k4[1] = h * dydx[1];

    y[0] += (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6;
    y[1] += (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6;
    x += h;

    x_vec.push_back(x);
    y_vec.push_back(y);
  }
}

void shooting_method(double x0, double xn, double h, double tol,
                     vector<double>& x_vec, vector<double>& y_vec) {
  double s0 = 1.0;
  double s1 = 3.0;

  auto F = [&](double s) -> double {
    vector<double> x_temp;
    vector<vector<double>> y_temp;
    vector<double> y0 = {s, 4.0};

    rk4_system(x0, y0, xn, h, x_temp, y_temp);

    double yN = y_temp.back()[0];   // y(1)
    double dyN = y_temp.back()[1];  // y'(1)
    return dyN - 2 * yN + 4;        // y'(1) - 2y(1) + 4 = 0
  };

  double F0 = F(s0);
  double F1 = F(s1);
  double s = s1;

  while (fabs(F1) > tol) {
    s = s1 - F1 * (s1 - s0) / (F1 - F0);
    s0 = s1;
    s1 = s;
    F0 = F1;
    F1 = F(s);
  }

  vector<vector<double>> y_temp_vec;
  rk4_system(x0, {s, 4.0}, xn, h, x_vec, y_temp_vec);

  y_vec.resize(x_vec.size());
  for (size_t i = 0; i < x_vec.size(); ++i) {
    y_vec[i] = y_temp_vec[i][0];
  }
}

void finite_difference(double x0, double xn, double h, vector<double>& x_vec,
                       vector<double>& y_vec) {
  int n = static_cast<int>((xn - x0) / h);
  x_vec.resize(n + 1);
  for (int i = 0; i <= n; ++i) {
    x_vec[i] = x0 + i * h;
  }

  vector<vector<double>> A(n + 1, vector<double>(n + 1, 0.0));
  vector<double> b(n + 1, 0.0);

  A[0][0] = -3.0;
  A[0][1] = 4.0;
  A[0][2] = -1.0;
  b[0] = 8 * h;

  A[n][n - 2] = 1.0;
  A[n][n - 1] = -4.0;
  A[n][n] = 3.0 - 4 * h;
  b[n] = -8 * h;

  for (int i = 1; i < n; ++i) {
    double x = x_vec[i];

    double alpha = 2 * x + h * (2 * x + 1);
    double beta = -4 * x + 4 * h * h;
    double gamma = 2 * x - h * (2 * x + 1);

    A[i][i - 1] = alpha;
    A[i][i] = beta;
    A[i][i + 1] = gamma;
    b[i] = 0.0;
  }

  y_vec.resize(n + 1, 0.0);

  for (int i = 0; i <= n; ++i) {
    double diag = A[i][i];
    for (int j = i; j <= n; ++j) {
      if (j > i)
        A[i][j] /= diag;
      else
        diag = A[i][i];
    }
    b[i] /= diag;

    for (int k = i + 1; k <= n; ++k) {
      double factor = A[k][i];
      for (int j = i; j <= n; ++j) {
        A[k][j] -= factor * A[i][j];
      }
      b[k] -= factor * b[i];
    }
  }

  y_vec[n] = b[n];
  for (int i = n - 1; i >= 0; --i) {
    y_vec[i] = b[i];
    for (int j = i + 1; j <= n; ++j) {
      y_vec[i] -= A[i][j] * y_vec[j];
    }
  }
}

vector<double> runge_romberg(const vector<double>& y_h,
                             const vector<double>& y_h2, int p) {
  vector<double> errors;
  for (size_t i = 0; i < y_h.size(); ++i) {
    size_t j = 2 * i;
    if (j < y_h2.size()) {
      double error = fabs(y_h[i] - y_h2[j]) / (pow(2, p) - 1);
      errors.push_back(error);
    }
  }
  return errors;
}

int main() {
  double x0 = 0.0;
  double xn = 1.0;
  double h = 0.1;
  double tol = 1e-6;

  vector<double> x_shoot, y_shoot;
  shooting_method(x0, xn, h, tol, x_shoot, y_shoot);

  vector<double> x_shoot_h2, y_shoot_h2;
  shooting_method(x0, xn, h / 2, tol, x_shoot_h2, y_shoot_h2);

  vector<double> rr_errors_shoot = runge_romberg(y_shoot, y_shoot_h2, 4);

  vector<double> x_fd, y_fd;
  finite_difference(x0, xn, h, x_fd, y_fd);

  vector<double> x_fd_h2, y_fd_h2;
  finite_difference(x0, xn, h / 2, x_fd_h2, y_fd_h2);

  vector<double> rr_errors_fd = runge_romberg(y_fd, y_fd_h2, 2);

  cout << "Shooting Method (h = 0.1)\n";
  cout << "x\t\tNumerical y\tExact y\t\tTrue Error\tRunge-Romberg Error\n";
  for (size_t i = 0; i < x_shoot.size(); ++i) {
    double exact = exact_y(x_shoot[i]);
    double true_error = fabs(y_shoot[i] - exact);
    double rr_error = (i < rr_errors_shoot.size()) ? rr_errors_shoot[i] : 0.0;

    cout << fixed << setprecision(1) << x_shoot[i] << fixed << setprecision(6)
         << "\t" << y_shoot[i] << "\t" << exact << "\t" << true_error << "\t"
         << rr_error << endl;
  }

  cout << "\nFinite Difference Method (h = 0.1)\n";
  cout << "x\t\tNumerical y\tExact y\t\tTrue Error\tRunge-Romberg Error\n";
  for (size_t i = 0; i < x_fd.size(); ++i) {
    double exact = exact_y(x_fd[i]);
    double true_error = fabs(y_fd[i] - exact);
    double rr_error = (i < rr_errors_fd.size()) ? rr_errors_fd[i] : 0.0;

    cout << fixed << setprecision(1) << x_fd[i] << fixed << setprecision(6)
         << "\t" << y_fd[i] << "\t" << exact << "\t" << true_error << "\t"
         << rr_error << endl;
  }

  return 0;
}