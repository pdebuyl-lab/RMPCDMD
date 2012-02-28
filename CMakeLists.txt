cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
project(MPCDMD)
FIND_PACKAGE(Git)
enable_language(Fortran)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

FIND_PACKAGE(HDF5_Fortran)

ENABLE_TESTING()

add_subdirectory(ParseText)
add_subdirectory(mtprng)

add_custom_target(MPCDMD_version ${CMAKE_COMMAND} -D
SRC=${CMAKE_CURRENT_SOURCE_DIR} -D
F95_MOD_DIR=${CMAKE_Fortran_MODULE_DIRECTORY} -D
CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH} -D
CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -D
CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS} -P
${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.cmake)

include_directories(${CMAKE_Fortran_MODULE_DIRECTORY} ${HDF5_INC}) 
add_executable(MPCDMD_test src/MPCDMD.f90 src/LJ.f90 src/MD.f90
src/MPCD.f90 src/group.f90 src/sys.f90 src/reaction.f90
src/polar_fields.f90
test/MPCDMD_test.f90 f90h5md/h5md.f90)
target_link_libraries(MPCDMD_test ParseText mtprng ${HDF5FORTRANLIB} ${HDF5LIB} ${ZLIB} ${SZLIB})
add_dependencies(MPCDMD_test MPCDMD_version)

add_executable(MPCD_test src/MPCDMD.f90 src/MPCD.f90 src/group.f90
src/sys.f90 test/MPCD_test.f90)
target_link_libraries(MPCD_test ParseText mtprng ${HDF5FORTRANLIB} ${HDF5LIB} ${ZLIB} ${SZLIB})
add_dependencies(MPCD_test MPCDMD_version)


