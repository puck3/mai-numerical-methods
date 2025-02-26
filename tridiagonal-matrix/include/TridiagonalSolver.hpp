#pragma once

#include <vector>

class TridiagonalSolver {
 private:
  std::vector<double> upperDiagonal;
  std::vector<double> mainDiagonal;
  std::vector<double> lowerDiagonal;
  std::vector<double> b;
  std::vector<double> x;
  int n;

  void validate();

 public:
  TridiagonalSolver(const std::vector<double>& lower,
                    const std::vector<double>& main,
                    const std::vector<double>& upper,
                    const std::vector<double>& rhs);

  void solve();

  const std::vector<double>& getSolution() const;
};