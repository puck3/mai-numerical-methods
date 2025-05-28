#include <cmath>
#include <iostream>
#include <vector>

double f(double x) { return sqrt(x) / (4 + 3 * x); }

double F(double x) {
  return ((2 * sqrt(x)) / 3) - ((4 * atan(sqrt(3 * x) / 2)) / (3 * sqrt(3)));
}

double rectangleMethod(double a, double b, double step) {
  double sum = 0;
  for (double x = a; x < b; x += step) {
    sum += f((x + x + step) / 2);
  }
  sum *= step;
  return sum;
}

double trapezoidalMethod(double a, double b, double step) {
  double sum = 0;
  for (double x = a; x < b; x += step) {
    sum += f(x) + f(x + step);
  }
  sum *= step / 2;
  return sum;
}

double simpsonMethod(double a, double b, double step) {
  double sum = 0;
  for (double x = a; x < b; x += step) {
    sum += f(x) + 4 * f(x + step / 2) + f(x + step);
  }
  sum *= step / 6;
  return sum;
}

double rungeRombergRichardsonMethod(double Fh, double Fkh, double k, double p) {
  return Fh + (Fh - Fkh) / (pow(k, p) - 1);
}

int main() {
  double a, b;
  std::cin >> a >> b;
  double Fabs = F(b) - F(a);

  double kh;
  std::cin >> kh;

  double FkhRectangle = rectangleMethod(a, b, kh);
  double FkhTrapezoidal = trapezoidalMethod(a, b, kh);
  double FkhSimpson = simpsonMethod(a, b, kh);

  std::cout << "Step = " << kh << std::endl;
  std::cout << "Rectangle method: " << FkhRectangle << std::endl;
  std::cout << "Trapezoidal method: " << FkhTrapezoidal << std::endl;
  std::cout << "Simpson method: " << FkhSimpson << std::endl;
  std::cout << std::endl;

  double h;
  std::cin >> h;

  double FhRectangle = rectangleMethod(a, b, h);
  double FhTrapezoidal = trapezoidalMethod(a, b, h);
  double FhSimpson = simpsonMethod(a, b, h);

  std::cout << "Step = " << h << std::endl;
  std::cout << "Rectangle method: " << FhRectangle << std::endl;
  std::cout << "Trapezoidal method: " << FhTrapezoidal << std::endl;
  std::cout << "Simpson method: " << FhSimpson << std::endl;
  std::cout << std::endl;

  std::cout << "Runge-Romberg Richardson Method:" << std::endl;

  double rectangleAccurate =
      rungeRombergRichardsonMethod(FhRectangle, FkhRectangle, kh / h, 2);
  std::cout << "Rectangle method: " << rectangleAccurate << std::endl;
  std::cout << "Error: " << fabs(Fabs - rectangleAccurate) << std::endl;

  double trapezoidalAccurate =
      rungeRombergRichardsonMethod(FhTrapezoidal, FkhTrapezoidal, kh / h, 2);
  std::cout << "Trapezoidal method: " << trapezoidalAccurate << std::endl;
  std::cout << "Error: " << fabs(Fabs - trapezoidalAccurate) << std::endl;

  double simpsonAccurate =
      rungeRombergRichardsonMethod(FhSimpson, FkhSimpson, kh / h, 4);
  std::cout << "Simpson method: " << simpsonAccurate << std::endl;
  std::cout << "Error: " << fabs(Fabs - simpsonAccurate) << std::endl;
}