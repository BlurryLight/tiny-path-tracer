cmake_minimum_required(VERSION 3.5)

project(Path_tracer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(OpenMP REQUIRED)
include_directories(${OpenMP_CXX_INCLUDE_DIRS})
add_executable(Path_tracer main.cpp vec3.h ray.h ray.cc hitable.h sphere.h hitableList.h camera.h)
target_link_libraries(Path_tracer PRIVATE OpenMP::OpenMP_CXX pthread)
