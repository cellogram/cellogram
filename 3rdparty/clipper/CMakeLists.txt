cmake_minimum_required(VERSION 3.1)
project(polyclipping)

add_library(polyclipping clipper/clipper.cpp)
target_include_directories(polyclipping PUBLIC .)
