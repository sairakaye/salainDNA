cmake_minimum_required(VERSION 3.15)
project(multicore_rm)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

add_executable(multicore_rm main.cpp common.cpp common.h pigeonhole.cpp pigeonhole.h edlib.cpp edlib.h searching.cpp searching.h minimizer.cpp minimizer.h bitmatrix.cpp bitmatrix.h)