#include <iostream>
#include <vector>

size_t find_left_index(const std::vector<double>& x, double value) {
  for (size_t i = 0; i < x.size(); ++i) {
    if (x[i] >= value) {
      return i - 1;
    }
  }
  return x.size() - 1;
}

double diff(double t, const std::vector<double>& x,
            const std::vector<double>& y) {
  size_t i = find_left_index(x, t);
  double res = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) +
               (((y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]) -
                 (y[i + 1] - y[i]) / (x[i + 1] - x[i])) /
                (x[i + 2] - x[i])) *
                   (t - x[i]) * (t - x[i + 1]);
  return res;
}

double diff2(double t, const std::vector<double>& x,
             const std::vector<double>& y) {
  size_t i = find_left_index(x, t);
  double res = 2 * (((y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]) -
                     (y[i + 1] - y[i]) / (x[i + 1] - x[i])) /
                    (x[i + 2] - x[i]));
  return res;
}

int main() {
  int n;
  std::cin >> n;
  std::vector<double> x(n);
  std::vector<double> y(n);
  for (size_t i = 0; i < n; ++i) {
    std::cin >> x[i];
  }
  for (size_t i = 0; i < n; ++i) {
    std::cin >> y[i];
  }
  double t;
  std::cin >> t;

  double y_1 = diff(t, x, y);
  double y_2 = diff2(t, x, y);
  std::cout << "y'(" << t << ") = " << y_1 << std::endl;
  std::cout << "y''(" << t << ") = " << y_2 << std::endl;
}