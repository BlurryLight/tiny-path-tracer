cmake_minimum_required(VERSION 3.5)

project(Path_tracer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(Threads REQUIRED)
file(COPY "config.ini" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set(HEADER_LIST
    headers/vec3.h headers/ray.h headers/hitable.h headers/sphere.h headers/hitable_list.h headers/camera.h
    headers/utils.h headers/aabb.h headers/texture.h headers/rect_box.h headers/material.h headers/perlin_noise.h
    third_party/inipp.h third_party/stb_image.h)

set(SRC_LIST
    src/ray.cc src/hitable.cc src/camera.cc src/hitable_list.cc src/sphere.cc src/utils.cc src/aabb.cc
    src/texture.cc src/rect_box.cc src/material.cc src/perlin_noise.cc
    )

IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    message("debug mode")
    add_definitions(-DDEBUG_MODE)
ENDIF(CMAKE_BUILD_TYPE MATCHES "Debug")

include_directories(headers/ third_party/)
add_executable(Path_tracer main.cpp ${HEADER_LIST} ${SRC_LIST})
target_link_libraries(Path_tracer PRIVATE Threads::Threads)
