cmake_minimum_required(VERSION 3.23)
project(omp-3)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

set(CMAKE_CXX_STANDARD 14)

add_executable(omp-3 main.cpp)
