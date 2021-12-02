cmake_minimum_required(VERSION 3.9)
# Get project name from folder name
get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
string(REPLACE " " "_" ProjectId ${ProjectId})
project(${ProjectId})
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")
# all .cc/.cpp files
set(SOURCES main.cc alpha_safe_paths.cc optimal_paths.cc draw_subgraph.cc)
add_executable(${ProjectId} ${SOURCES})
# compiler flags for debugging
if (APPLE)
    set(CMAKE_CXX_FLAGS "-O3 -openmp")
else()
    set(CMAKE_CXX_FLAGS "-O3 -fopenmp")
endif()

# GMP
find_package(GMP REQUIRED)
# if .h files are in other folder?
include_directories(./ ${GMP_INCLUDE_DIR})
target_link_libraries(${ProjectId} ${GMP_LIBRARY})
# OpenMP
#find_package(OpenMP REQUIRED)
#target_link_libraries(${ProjectId} OpenMP::OpenMP_CXX)
# runnable output file
set_target_properties(${ProjectId} PROPERTIES OUTPUT_NAME "main")