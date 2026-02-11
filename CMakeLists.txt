cmake_minimum_required(VERSION 3.16)
project(your_project LANGUAGES CXX)

find_package(Eigen3 REQUIRED)

add_library(core INTERFACE)
target_include_directories(core INTERFACE include src)
target_link_libraries(core INTERFACE Eigen3::Eigen)

add_subdirectory(tests)
