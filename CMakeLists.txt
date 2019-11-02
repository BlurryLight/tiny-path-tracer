cmake_minimum_required(VERSION 3.5)

project(Path_tracer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(Path_tracer main.cpp vec3.h)
