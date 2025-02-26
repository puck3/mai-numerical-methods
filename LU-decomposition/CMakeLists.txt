cmake_minimum_required(VERSION 3.28)
project(LU)

include_directories(include)

add_executable(LU src/main.cpp src/Matrix.cpp src/LUDecomposer.cpp)
