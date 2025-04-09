#include "Matrix.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

Matrix::Matrix(size_t n, size_t m, double init)
    : rows(n), cols(m), data(n, std::vector<double>(m, init)) {}

Matrix::Matrix(const std::vector<std::vector<double>>& d)
    : data(d), rows(d.size()), cols(d.empty() ? 0 : d[0].size()) {}

size_t Matrix::getRows() const { return rows; }
size_t Matrix::getCols() const { return cols; }

double& Matrix::operator()(size_t i, size_t j) { return data[i][j]; }

const double& Matrix::operator()(size_t i, size_t j) const {
  return data[i][j];
}

Matrix Matrix::transpose() const {
  Matrix result(cols, rows);
  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j) result(j, i) = data[i][j];
  return result;
}

Matrix Matrix::identity(size_t n) {
  Matrix result(n, n);
  for (size_t i = 0; i < n; ++i) result(i, i) = 1.0;
  return result;
}

void Matrix::print(const std::string& title) const {
  if (!title.empty()) std::cout << title << "\n";
  for (const auto& row : data) {
    for (double val : row)
      std::cout << std::setw(12) << std::setprecision(6) << val << " ";
    std::cout << "\n";
  }
}

Matrix Matrix::operator*(const Matrix& other) const {
  if (cols != other.rows)
    throw std::invalid_argument("Matrix dimensions mismatch");

  Matrix result(rows, other.cols);
  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < other.cols; ++j)
      for (size_t k = 0; k < cols; ++k)
        result(i, j) += data[i][k] * other(k, j);
  return result;
}