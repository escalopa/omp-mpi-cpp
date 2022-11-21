cmake_minimum_required(VERSION 3.23)
project(mpi-0)

set(CMAKE_CXX_STANDARD 14)

include_directories(/usr/include/mpi/)

set(CMAKE_C_COMPILER "/usr/bin/mpicc")
set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")

add_executable(mpi-0 main.c)
