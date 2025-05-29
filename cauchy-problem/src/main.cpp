#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

double exact_y(double x) { return pow(x, 4) + pow(x, -3); }

double exact_dy(double x) { return 4 * pow(x, 3) - 3 * pow(x, -4); }

double f_u(double x, double u, double v) { return v; }

double f_v(double x, double u, double v) { return 12.0 * u / (x * x); }

void euler(double x0, double u0, double v0, double xn, double h,
           vector<double>& x_vec, vector<double>& u_vec,
           vector<double>& v_vec) {
  x_vec.clear();
  u_vec.clear();
  v_vec.clear();
  double x = x0;
  double u = u0;
  double v = v0;

  x_vec.push_back(x);
  u_vec.push_back(u);
  v_vec.push_back(v);

  while (x < xn) {
    double u_next = u + h * f_u(x, u, v);
    double v_next = v + h * f_v(x, u, v);
    x += h;

    u = u_next;
    v = v_next;

    x_vec.push_back(x);
    u_vec.push_back(u);
    v_vec.push_back(v);
  }
}

void rk4(double x0, double u0, double v0, double xn, double h,
         vector<double>& x_vec, vector<double>& u_vec, vector<double>& v_vec) {
  x_vec.clear();
  u_vec.clear();
  v_vec.clear();
  double x = x0;
  double u = u0;
  double v = v0;

  x_vec.push_back(x);
  u_vec.push_back(u);
  v_vec.push_back(v);

  while (x < xn) {
    double k1u = h * f_u(x, u, v);
    double k1v = h * f_v(x, u, v);

    double k2u = h * f_u(x + h / 2, u + k1u / 2, v + k1v / 2);
    double k2v = h * f_v(x + h / 2, u + k1u / 2, v + k1v / 2);

    double k3u = h * f_u(x + h / 2, u + k2u / 2, v + k2v / 2);
    double k3v = h * f_v(x + h / 2, u + k2u / 2, v + k2v / 2);

    double k4u = h * f_u(x + h, u + k3u, v + k3v);
    double k4v = h * f_v(x + h, u + k3u, v + k3v);

    u += (k1u + 2 * k2u + 2 * k3u + k4u) / 6;
    v += (k1v + 2 * k2v + 2 * k3v + k4v) / 6;
    x += h;

    x_vec.push_back(x);
    u_vec.push_back(u);
    v_vec.push_back(v);
  }
}

void adams4(double x0, double u0, double v0, double xn, double h,
            vector<double>& x_vec, vector<double>& u_vec,
            vector<double>& v_vec) {
  x_vec.clear();
  u_vec.clear();
  v_vec.clear();
  vector<double> f_u_vals, f_v_vals;

  vector<double> x_temp, u_temp, v_temp;
  rk4(x0, u0, v0, x0 + 3 * h, h, x_temp, u_temp, v_temp);

  for (int i = 0; i < 4; ++i) {
    x_vec.push_back(x_temp[i]);
    u_vec.push_back(u_temp[i]);
    v_vec.push_back(v_temp[i]);
    f_u_vals.push_back(f_u(x_temp[i], u_temp[i], v_temp[i]));
    f_v_vals.push_back(f_v(x_temp[i], u_temp[i], v_temp[i]));
  }

  double x = x0 + 3 * h;
  int idx = 3;

  while (x < xn - 1e-9) {
    double fu_next =
        f_u_vals[idx] + h / 24 *
                            (55 * f_u_vals[idx] - 59 * f_u_vals[idx - 1] +
                             37 * f_u_vals[idx - 2] - 9 * f_u_vals[idx - 3]);
    double fv_next =
        f_v_vals[idx] + h / 24 *
                            (55 * f_v_vals[idx] - 59 * f_v_vals[idx - 1] +
                             37 * f_v_vals[idx - 2] - 9 * f_v_vals[idx - 3]);

    double u_next = u_vec[idx] + h * fu_next;
    double v_next = v_vec[idx] + h * fv_next;
    x += h;
    idx++;

    x_vec.push_back(x);
    u_vec.push_back(u_next);
    v_vec.push_back(v_next);

    f_u_vals.push_back(f_u(x, u_next, v_next));
    f_v_vals.push_back(f_v(x, u_next, v_next));
  }
}

void runge_romberg(vector<double>& u_h, vector<double>& u_h2, int p,
                   vector<double>& errors) {
  errors.resize(u_h.size());
  for (size_t i = 0; i < u_h.size(); ++i) {
    errors[i] = abs(u_h[i] - u_h2[2 * i]) / (pow(2, p) - 1);
  }
}

int main() {
  double x0 = 1.0;
  double u0 = 2.0;  // y(1)
  double v0 = 1.0;  // y'(1)
  double xn = 2.0;
  double h = 0.1;
  double h2 = h / 2;

  vector<double> x_euler, u_euler, v_euler;
  euler(x0, u0, v0, xn, h, x_euler, u_euler, v_euler);

  vector<double> x_euler_h2, u_euler_h2, v_euler_h2;
  euler(x0, u0, v0, xn, h2, x_euler_h2, u_euler_h2, v_euler_h2);

  vector<double> rr_errors_euler;
  runge_romberg(u_euler, u_euler_h2, 1, rr_errors_euler);

  vector<double> x_rk4, u_rk4, v_rk4;
  rk4(x0, u0, v0, xn, h, x_rk4, u_rk4, v_rk4);

  vector<double> x_rk4_h2, u_rk4_h2, v_rk4_h2;
  rk4(x0, u0, v0, xn, h2, x_rk4_h2, u_rk4_h2, v_rk4_h2);

  vector<double> rr_errors_rk4;
  runge_romberg(u_rk4, u_rk4_h2, 4, rr_errors_rk4);

  vector<double> x_adams, u_adams, v_adams;
  adams4(x0, u0, v0, xn, h, x_adams, u_adams, v_adams);

  vector<double> x_adams_h2, u_adams_h2, v_adams_h2;
  adams4(x0, u0, v0, xn, h2, x_adams_h2, u_adams_h2, v_adams_h2);

  vector<double> rr_errors_adams;
  runge_romberg(u_adams, u_adams_h2, 4, rr_errors_adams);

  cout << "Euler Method (h = 0.1)\n";
  cout << "x\t\tNumerical y\tExact y\t\tTrue Error\tRunge-Romberg Error\n";
  for (size_t i = 0; i < x_euler.size(); ++i) {
    double exact = exact_y(x_euler[i]);
    double true_error = abs(u_euler[i] - exact);
    cout << fixed << setprecision(1) << x_euler[i] << fixed << setprecision(6)
         << "\t" << u_euler[i] << "\t" << exact << "\t" << true_error << "\t"
         << rr_errors_euler[i] << endl;
  }

  cout << "\nRunge-Kutta 4th Order (h = 0.1)\n";
  cout << "x\t\tNumerical y\tExact y\t\tTrue Error\tRunge-Romberg Error\n";
  for (size_t i = 0; i < x_rk4.size(); ++i) {
    double exact = exact_y(x_rk4[i]);
    double true_error = abs(u_rk4[i] - exact);
    cout << fixed << setprecision(1) << x_rk4[i] << fixed << setprecision(6)
         << "\t" << u_rk4[i] << "\t" << exact << "\t" << true_error << "\t"
         << rr_errors_rk4[i] << endl;
  }

  cout << "\nAdams 4th Order (h = 0.1)\n";
  cout << "x\t\tNumerical y\tExact y\t\tTrue Error\tRunge-Romberg Error\n";
  for (size_t i = 0; i < x_adams.size(); ++i) {
    double exact = exact_y(x_adams[i]);
    double true_error = abs(u_adams[i] - exact);
    cout << fixed << setprecision(1) << x_adams[i] << fixed << setprecision(6)
         << "\t" << u_adams[i] << "\t" << exact << "\t" << true_error << "\t"
         << rr_errors_adams[i] << endl;
  }

  return 0;
}