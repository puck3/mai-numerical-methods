#include <iostream>
#include <stdexcept>
#include <vector>

#include "LUDecomposer.hpp"
#include "Matrix.hpp"

using namespace std;

int main() {
  try {
    size_t n;
    cin >> n;
    vector<vector<double>> MatrixData(n, vector<double>(n));
    vector<double> b(n);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        cin >> MatrixData[i][j];
      }
      cin >> b[i];
    }
    Matrix A(MatrixData);
    LUDecomposer lu(A);

    vector<double> x = lu.solve(b);
    cout << "Solution x: ";
    for (double xi : x) cout << xi << " ";
    cout << endl;

    cout << "Determinant: " << lu.determinant() << endl;

    Matrix inv = lu.inverse();
    inv.print("Inverse matrix:");

  } catch (const exception& e) {
    cerr << "Error: " << e.what() << endl;
  }
  return 0;
}