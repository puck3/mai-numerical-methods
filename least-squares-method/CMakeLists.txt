cmake_minimum_required(VERSION 3.28)
project(LSM)

include_directories(include)

add_executable(LSM src/main.cpp src/Matrix.cpp src/LUDecomposer.cpp)
