cmake_minimum_required(VERSION 3.5)

project(Path_tracer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(Threads REQUIRED)
file(COPY "config.ini" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set(HEADER_LIST
    headers/vec3.h headers/ray.h headers/hitable.h headers/sphere.h headers/hitable_list.h headers/camera.h
    headers/utils.h headers/aabb.h headers/texture.h
    third_party/inipp.h)

set(SRC_LIST
    src/ray.cc src/hitable.cc src/camera.cc src/hitable_list.cc src/sphere.cc src/utils.cc src/aabb.cc
    src/texture.cc
    )
include_directories(headers/)
add_executable(Path_tracer main.cpp ${HEADER_LIST} ${SRC_LIST})
target_link_libraries(Path_tracer PRIVATE Threads::Threads)
