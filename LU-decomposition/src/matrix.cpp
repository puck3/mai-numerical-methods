#include "matrix.hpp"

Matrix::Matrix(size_t n, size_t m)
    : rows(n),
      cols(m),
      data(std::vector<std::vector<double>>(rows,
                                            std::vector<double>(cols, 0))) {}

double& Matrix::operator()(size_t i, size_t j) { return this->data[i][j]; }

const double& Matrix::operator()(size_t i, size_t j) const {
  return this->data[i][j];
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) {
  for (size_t i = 0; i < m.rows; ++i) {
    for (size_t j = 0; j < m.cols; ++j) {
      os << m(i, j) << " ";
    }
    os << std::endl;
  }
  return os;
}

std::istream& operator>>(std::istream& is, Matrix& m) {
  for (size_t i = 0; i < m.rows; ++i) {
    for (size_t j = 0; j < m.cols; ++j) {
      is >> m(i, j);
    }
  }
  return is;
}