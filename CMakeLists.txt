cmake_minimum_required(VERSION 3.28)
project(numerical-methods)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(LU-decomposition)
add_subdirectory(tridiagonal-matrix)
add_subdirectory(iterative-methods)
