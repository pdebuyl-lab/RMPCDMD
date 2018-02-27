cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_SOURCE_DIR}/cmake/compiler.cmake)
project(RMPCDMD Fortran)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

enable_testing()

set(PT_HDF5 ON CACHE BOOL "")
add_subdirectory(ParseText)
add_subdirectory(fortran_h5md)
add_subdirectory(fortran_tester)
add_subdirectory(random_module)
add_subdirectory(fortran_quaternion)

find_package(Git)
find_package(HDF5 COMPONENTS Fortran REQUIRED)

include_directories(${CMAKE_Fortran_MODULE_DIRECTORY}
${HDF5_Fortran_INCLUDE_DIR} ${HDF5_Fortran_INCLUDE_DIRS})

set(RMPCDMD_PROGRAMS poiseuille_flow single_dimer_pbc single_dimer single_body
chemotactic_cell single_sphere_thermo_trap
n_colloids_pbc)

configure_file(scripts/rmpcdmd.in rmpcdmd @ONLY)
install(PROGRAMS ${CMAKE_BINARY_DIR}/rmpcdmd DESTINATION bin)

add_custom_target(RMPCDMD_REVISION_TARGET ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.cmake)
execute_process(COMMAND ${CMAKE_COMMAND} -D RMPCDMD_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.cmake)

add_library(RMPCDMD
src/common.f90
src/correlator.f90
src/cell_system.f90
src/hilbert.f90
src/interaction.f90
src/particle_system.f90
src/particle_system_io.f90
src/mpcd.f90
src/md.f90
src/neighbor_list.f90
src/polar_fields.f90
src/planar_fields.f90
${RMPCDMD_BINARY_DIR}/rmpcdmd_module.f90)
target_link_libraries(RMPCDMD fortran_h5md mt19937ar_module threefry_module quaternion)
add_dependencies(RMPCDMD RMPCDMD_REVISION_TARGET)

add_executable(setup_simple_fluid
test/setup_simple_fluid.f90)
target_link_libraries(setup_simple_fluid RMPCDMD fortran_h5md ParseText)

add_executable(setup_bulk_decay
test/setup_bulk_decay.f90)
target_link_libraries(setup_bulk_decay RMPCDMD fortran_h5md ParseText)

add_executable(setup_thermal_fluid
test/setup_thermal_fluid.f90)
target_link_libraries(setup_thermal_fluid RMPCDMD fortran_h5md ParseText)

foreach(program_name IN LISTS RMPCDMD_PROGRAMS)
add_executable(${program_name} programs/${program_name}.f90)
target_link_libraries(${program_name} RMPCDMD fortran_h5md ParseText)
endforeach(program_name)

add_executable(setup_simple_colloids
test/setup_simple_colloids.f90)
target_link_libraries(setup_simple_colloids RMPCDMD fortran_h5md ParseText)

add_executable(setup_single_colloid
test/setup_single_colloid.f90)
target_link_libraries(setup_single_colloid RMPCDMD fortran_h5md ParseText)

add_executable(setup_simple_rattle
test/setup_simple_rattle.f90)
target_link_libraries(setup_simple_rattle RMPCDMD fortran_h5md ParseText)

add_executable(setup_single_catalytic_fixed_sphere
test/setup_single_catalytic_fixed_sphere.f90)
target_link_libraries(setup_single_catalytic_fixed_sphere RMPCDMD fortran_h5md ParseText)

add_executable(test_cell_system_0
test/test_cell_system_0.f90)
target_link_libraries(test_cell_system_0 RMPCDMD fortran_tester)

add_executable(test_common_0
test/test_common_0.f90)
target_link_libraries(test_common_0 RMPCDMD fortran_tester)

add_executable(test_common_1
test/test_common_1.f90)
target_link_libraries(test_common_1 RMPCDMD fortran_tester)

add_executable(test_neighbor_list_0
test/test_neighbor_list_0.f90)
target_link_libraries(test_neighbor_list_0 RMPCDMD fortran_tester)

add_executable(test_neighbor_list_1
test/test_neighbor_list_1.f90)
target_link_libraries(test_neighbor_list_1 RMPCDMD fortran_tester)

add_executable(test_particle_system_0
test/test_particle_system_0.f90)
target_link_libraries(test_particle_system_0 RMPCDMD fortran_tester)

add_executable(test_correlator_0
test/test_correlator_0.f90)
target_link_libraries(test_correlator_0 RMPCDMD fortran_tester)

add_executable(test_elastic_network_0
test/test_elastic_network_0.f90)
target_link_libraries(test_elastic_network_0 RMPCDMD fortran_tester)

add_executable(test_md_0
test/test_md_0.f90)
target_link_libraries(test_md_0 RMPCDMD fortran_tester)

add_test(NAME test_cell_system_0 COMMAND test_cell_system_0)
add_test(NAME test_common_0 COMMAND test_common_0)
add_test(NAME test_common_1 COMMAND test_common_1)
add_test(NAME test_neighbor_list_0 COMMAND test_neighbor_list_0)
add_test(NAME test_neighbor_list_1 COMMAND test_neighbor_list_1)
add_test(NAME test_particle_system_0 COMMAND test_particle_system_0)
add_test(NAME test_correlator_0 COMMAND test_correlator_0)
add_test(NAME test_elastic_network_0 COMMAND test_elastic_network_0)
add_test(NAME test_md_0 COMMAND test_md_0)

set_tests_properties(test_neighbor_list_1 PROPERTIES WILL_FAIL TRUE)
