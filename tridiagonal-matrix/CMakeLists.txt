cmake_minimum_required(VERSION 3.28)
project(tridiagonal)

include_directories(include)

add_executable(tridiagonal src/main.cpp src/TridiagonalSolver.cpp)
