#pragma once

#include <vector>

class TridiagonalSolver {
 private:
  std::vector<double> a;
  std::vector<double> b;
  std::vector<double> c;
  std::vector<double> d;
  std::vector<double> x;
  int n;

  void validate();

 public:
  TridiagonalSolver(const std::vector<double>& a, const std::vector<double>& b,
                    const std::vector<double>& c, const std::vector<double>& d);

  void solve();

  const std::vector<double>& getSolution() const;
};