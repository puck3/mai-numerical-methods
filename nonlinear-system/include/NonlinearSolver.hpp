
#include <cmath>
#include <stdexcept>
#include <vector>

struct Result {
  std::vector<double> errors;
  double x1;
  double x2;
  int iterations;
};

class NonlinearSolver {
 private:
  double epsilon;
  int maxIterations;

  double f1(double x1, double x2) const { return x1 * x1 - 2 * log10(x2) - 1; }

  double f2(double x1, double x2) const { return x1 * x1 - 2 * x1 * x2 + 2; }

  void computeJacobian(double x1, double x2, double& df1_dx1, double& df1_dx2,
                       double& df2_dx1, double& df2_dx2) const {
    df1_dx1 = 2 * x1;
    df1_dx2 = -2 / (x2 * log(10.0));
    df2_dx1 = 2 * x1 - 2 * x2;
    df2_dx2 = -2 * x1;
  }

  double norm(double a, double b) const { return sqrt(a * a + b * b); }

 public:
  NonlinearSolver(double eps, int maxIter)
      : epsilon(eps), maxIterations(maxIter) {}

  Result simpleIteration(double x1, double x2) const {
    double q = 0.55;
    Result result;
    for (int i = 0; i < maxIterations; ++i) {
      double x1_new = sqrt(2 * log10(x2) + 1);
      double x2_new = (x1 * x1 + 2) / (2 * x1);

      double error = (q / (1 - q)) * norm(x1_new - x1, x2_new - x2);
      result.errors.push_back(error);

      double res = norm(f1(x1_new, x2_new), f2(x1_new, x2_new));

      x1 = x1_new;
      x2 = x2_new;

      if (error <= epsilon) break;
    }

    result.x1 = x1;
    result.x2 = x2;
    result.iterations = result.errors.size();
    return result;
  }

  Result newtonMethod(double x1_initial, double x2_initial) const {
    Result result;
    double x1 = x1_initial;
    double x2 = x2_initial;

    for (int i = 0; i < maxIterations; ++i) {
      double f1_val = f1(x1, x2);
      double f2_val = f2(x1, x2);

      double df1_dx1, df1_dx2, df2_dx1, df2_dx2;
      computeJacobian(x1, x2, df1_dx1, df1_dx2, df2_dx1, df2_dx2);

      double det = df1_dx1 * df2_dx2 - df1_dx2 * df2_dx1;
      if (det == 0) throw std::runtime_error("Jacobian matrix is singular");

      double delta1 = (-f1_val * df2_dx2 + f2_val * df1_dx2) / det;
      double delta2 = (f1_val * df2_dx1 - f2_val * df1_dx1) / det;

      double x1_new = x1 + delta1;
      double x2_new = x2 + delta2;

      double error = norm(delta1, delta2);
      result.errors.push_back(error);

      double res = norm(f1(x1_new, x2_new), f2(x1_new, x2_new));

      x1 = x1_new;
      x2 = x2_new;

      if (error <= epsilon) break;
    }

    result.x1 = x1;
    result.x2 = x2;
    result.iterations = result.errors.size();
    return result;
  }
};