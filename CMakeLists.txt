cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_SOURCE_DIR}/cmake/compiler.cmake)
project(RMPCDMD Fortran)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

enable_testing()

add_subdirectory(ParseText)
add_subdirectory(fortran_h5md)
add_subdirectory(fortran_tester)
add_subdirectory(random_module)

include_directories(${CMAKE_Fortran_MODULE_DIRECTORY} ${HDF5_Fortran_INCLUDE_DIR})

add_library(RMPCDMD
src/common.f90
src/cell_system.f90
src/hilbert.f90
src/interaction.f90
src/particle_system.f90
src/mpcd.f90
src/neighbor_list.f90)
target_link_libraries(RMPCDMD fortran_h5md mt19937ar_module)

add_executable(try_all
test/try_all.f90)
target_link_libraries(try_all RMPCDMD fortran_h5md ParseText)

add_executable(setup_simple_fluid
test/setup_simple_fluid.f90)
target_link_libraries(setup_simple_fluid RMPCDMD fortran_h5md ParseText)

add_executable(setup_simple_colloids
test/setup_simple_colloids.f90)
target_link_libraries(setup_simple_colloids RMPCDMD fortran_h5md ParseText)

add_executable(test_cell_system_0
test/test_cell_system_0.f90)
target_link_libraries(test_cell_system_0 RMPCDMD fortran_tester)

add_test(test_cell_system_0 test_cell_system_0)
set_tests_properties(test_cell_system_0 PROPERTIES PASS_REGULAR_EXPRESSION "all tests succeeded")
