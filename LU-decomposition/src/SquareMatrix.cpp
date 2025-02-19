#include "SquareMatrix.hpp"

SquareMatrix::SquareMatrix(size_t n)
    : n(n),
      data(std::vector<std::vector<double>>(n, std::vector<double>(n, 0))) {}

size_t SquareMatrix::size() const { return n; }

double& SquareMatrix::operator()(size_t i, size_t j) {
  return this->data[i][j];
}

const double& SquareMatrix::operator()(size_t i, size_t j) const {
  return this->data[i][j];
}

std::ostream& operator<<(std::ostream& os, const SquareMatrix& m) {
  for (size_t i = 0; i < m.size(); ++i) {
    for (size_t j = 0; j < m.size(); ++j) {
      os << m(i, j) << " ";
    }
    os << std::endl;
  }
  return os;
}

std::istream& operator>>(std::istream& is, SquareMatrix& m) {
  for (size_t i = 0; i < m.size(); ++i) {
    for (size_t j = 0; j < m.size(); ++j) {
      is >> m(i, j);
    }
  }
  return is;
}