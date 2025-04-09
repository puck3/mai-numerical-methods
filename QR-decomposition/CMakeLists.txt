cmake_minimum_required(VERSION 3.28)
project(QR)

include_directories(include)
add_executable(QR src/main.cpp src/Matrix.cpp src/QRDecomposer.cpp src/QREigenSolver.cpp)
