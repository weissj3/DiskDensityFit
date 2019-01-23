project(DiskFit)

cmake_minimum_required(VERSION 2.6)


set(CMAKE_CXX_FLAGS "-g -Wall -O3 -funroll-loops -ftree-vectorize -msse3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/local/lib)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/undvc_common/cmake/Modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules/")



find_package(Boost REQUIRED COMPONENTS filesystem system)


add_subdirectory(tao)
add_subdirectory(undvc_common)
add_subdirectory(Tests)
add_subdirectory(Density_Program)
include_directories(${PROJECT_SOURCE_DIR}/tao)
include_directories(${PROJECT_SOURCE_DIR}/undvc_common)
include_directories(${PROJECT_SOURCE_DIR}/Tests)
include_directories(${PROJECT_SOURCE_DIR}/Density_Program)


