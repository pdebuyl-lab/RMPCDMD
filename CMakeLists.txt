cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
project(RMPCDMD C Fortran)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

enable_testing()

add_subdirectory(ParseText)

include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
add_executable(try_all
src/common.f90
src/cell_system.f90
src/particle_system.f90
src/hilbert.f90
src/neighbor_list.f90
test/try_all.f90)
target_link_libraries(try_all ParseText dl)

add_executable(test_cell_system_0
src/cell_system.f90
src/hilbert.f90
test/test_cell_system_0.f90)
target_link_libraries(test_cell_system_0 dl)

add_test(test_cell_system_0 test_cell_system_0)
set_tests_properties(test_cell_system_0 PROPERTIES PASS_REGULAR_EXPRESSION "success")
