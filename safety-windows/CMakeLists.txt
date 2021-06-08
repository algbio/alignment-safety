cmake_minimum_required(VERSION 2.8)
# Get project name from folder name
get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
string(REPLACE " " "_" ProjectId ${ProjectId})
project(${ProjectId})
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")
find_package(GMP REQUIRED)
# all .cc/.cpp files
set(SOURCES main.cc safety_windows.cc alpha_safe_paths.cc optimal_paths.cc)
add_executable(${PROJECT_NAME} ${SOURCES})
# compiler flags for debugging
set(CMAKE_CXX_FLAGS "-g -W -O0")
# if .h files are in other folder?
include_directories(./ ${GMP_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} gmp libgmp.so)
# runnable output file
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "main")
