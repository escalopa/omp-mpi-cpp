cmake_minimum_required(VERSION 3.23)
project(omp-0)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")

add_executable(omp-0 main.cpp)
