cmake_minimum_required(VERSION 3.16)
project(project2 LANGUAGES CXX)

add_library(core INTERFACE)
target_include_directories(core INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    /sw/pkgs/arc/eigen/eigen-3.4.0
)

add_subdirectory(tests)
