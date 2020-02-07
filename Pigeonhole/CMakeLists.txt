cmake_minimum_required(VERSION 3.15)
project(Pigeonhole)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

add_executable(Pigeonhole main.cpp common.cpp common.h pigeonhole.h pigeonhole.cpp minimizers.cpp minimizers.h addressing.cpp addressing.h edlib.cpp edlib.h)