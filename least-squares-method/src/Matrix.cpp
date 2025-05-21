#include "Matrix.hpp"

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

void Matrix::swap_rows(size_t i, size_t j) { swap(data[i], data[j]); }

std::vector<double>& Matrix::operator[](size_t i) { return data[i]; }
const std::vector<double>& Matrix::operator[](size_t i) const {
  return data[i];
}

void Matrix::print(const std::string& title) const {
  if (!title.empty()) std::cout << title << std::endl;
  for (const auto& row : data) {
    for (double val : row) std::cout << val << "\t";
    std::cout << std::endl;
  }
}