#include <iostream>
#include <vector>

#include "matrix.hpp"

int main() {
  int n = 4;
  Matrix m(4, 4);
  std::cin >> m;
  std::cout << m;
}