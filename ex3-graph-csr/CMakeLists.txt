cmake_minimum_required(VERSION 3.17)
project(ex3_graph_csr CXX)

set(CMAKE_CXX_STANDARD 17)

add_executable(ex3_graph_csr "main.cpp")

find_package(OpenMP REQUIRED)
target_link_libraries(ex3_graph_csr PRIVATE OpenMP::OpenMP_CXX)