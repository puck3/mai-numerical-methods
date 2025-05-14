
#include <vector>

class LagrangePolynomial {
 private:
  std::vector<double> f;
  std::vector<double> x;

 public:
  LagrangePolynomial(double (*func)(double), const std::vector<double>& nodes)
      : f(), x(nodes) {
    for (size_t i = 0; i < x.size(); ++i) {
      f.push_back(func(x[i]));
    }
  }

  double interpolate(double t) {
    double res = 0;
    for (size_t i = 0; i < x.size(); ++i) {
      double prod = f[i];
      for (size_t j = 0; j < x.size(); ++j) {
        if (j == i) {
          continue;
        }
        prod *= (t - x[j]) / (x[i] - x[j]);
      }
      res += prod;
    }
    return res;
  }
};