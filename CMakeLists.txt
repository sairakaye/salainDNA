cmake_minimum_required(VERSION 3.10)
project(salainDNA)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

add_executable(salainDNA main.cpp common.cpp common.h seedselector.cpp seedselector.h edlib.cpp edlib.h searching.cpp searching.h minimizers.cpp minimizers.h bitmatrix.cpp bitmatrix.h command.cpp command.h indexing.cpp indexing.h directaddressing.cpp directaddressing.h openaddressing.cpp openaddressing.h output.cpp output.h searchall.cpp searchall.h searchexit.cpp searchexit.h verification.cpp verification.h)