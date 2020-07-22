cmake_minimum_required(VERSION 3.15)
project(multicore_rm)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

add_executable(multicore_rm main.cpp common.cpp common.h pigeonhole.cpp pigeonhole.h edlib.cpp edlib.h searching.cpp searching.h minimizers.cpp minimizers.h bitmatrix.cpp bitmatrix.h command.cpp command.h indexing.cpp indexing.h directaddressing.cpp directaddressing.h openaddressing.cpp openaddressing.h output.cpp output.h)