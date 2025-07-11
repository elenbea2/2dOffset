cmake_minimum_required(VERSION 3.10)
project(offset_project)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CGAL REQUIRED)

add_executable(offset offset_holes.cpp)

target_link_libraries(offset CGAL::CGAL)
