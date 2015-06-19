cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
project(RMPCDMD C Fortran)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

add_subdirectory(ParseText)

include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
add_executable(try_all
src/cell_system.f90
src/particle_system.f90
src/hilbert.f90
test/try_all.f90)
target_link_libraries(try_all ParseText dl)

