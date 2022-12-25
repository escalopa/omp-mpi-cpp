cmake_minimum_required(VERSION 3.24)
project(f2dtc)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_COMPILER "/usr/bin/mpicc")
set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")

add_executable(f2dtc main.cpp)
