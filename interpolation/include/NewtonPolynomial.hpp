#include <vector>

class NewtonPolynomial {
 private:
  std::vector<double> f;
  std::vector<double> x;

 public:
  NewtonPolynomial(double (*func)(double), const std::vector<double>& nodes)
      : f(), x(nodes) {
    for (size_t i = 0; i < x.size(); ++i) {
      f.push_back(func(x[i]));
    }

    for (size_t j = 1; j < x.size(); ++j) {
      for (size_t i = x.size() - 1; i >= j; --i) {
        f[i] = (f[i] - f[i - 1]) / (x[i] - x[i - j]);
      }
    }
  }

  double interpolate(double t) {
    double res = 0;
    double coef = 1;
    for (size_t i = 0; i < x.size(); ++i) {
      res += coef * f[i];
      coef *= (t - x[i]);
    }
    return res;
  }
};