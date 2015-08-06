set(CMAKE_BUILD_TYPE_INIT "Release")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_INIT "-fopenmp -Wall -Wextra -Wconversion -std=f2008 -pedantic")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS_INIT "-openmp -stand f08")
endif()
