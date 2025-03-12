#pragma once

#include <vector>

#include "IterativeSolver.hpp"

using namespace std;

class SeidelSolver : public IterativeSolver {
 public:
  using IterativeSolver::IterativeSolver;

  vector<double> solve() override {
    prepare();
    int n = A.size();
    vector<double> x_prev(n, 0.0);
    iterations = 0;

    do {
      x_prev = x;
      for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
          if (i != j) sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
      }
      ++iterations;

      vector<double> diff(n);
      for (int i = 0; i < n; ++i) diff[i] = x[i] - x_prev[i];
      if (norm(diff) < epsilon) break;

    } while (iterations < max_iterations);

    return x;
  }
};