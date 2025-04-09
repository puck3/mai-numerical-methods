#pragma once

#include <vector>
#include <iostream>

using namespace std;

class Matrix {
 private:
  vector<vector<double>> data;
  size_t rows, cols;

 public:
  Matrix(size_t n = 0, size_t m = 0, double init = 0.0)
      : rows(n), cols(m), data(n, vector<double>(m, init)) {}

  Matrix(const vector<vector<double>>& d)
      : data(d), rows(d.size()), cols(d.empty() ? 0 : d[0].size()) {}

  size_t getRows() const { return rows; }
  size_t getCols() const { return cols; }

  double& operator()(size_t i, size_t j) { return data[i][j]; }
  const double& operator()(size_t i, size_t j) const { return data[i][j]; }

  void print(const string& title = "") const {
    if (!title.empty()) cout << title << endl;
    for (const auto& row : data) {
      for (double val : row) printf("%10.6f ", val);
      cout << endl;
    }
  }
};