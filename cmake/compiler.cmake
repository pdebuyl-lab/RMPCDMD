set(CMAKE_BUILD_TYPE_INIT "Release")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_INIT "-fopenmp -Wall -Wextra -Wconversion -std=f2008 -pedantic")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS_INIT "-openmp -stand f08")
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS_INIT "-Wall -std=c99")
elseif(CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_C_FLAGS_INIT "-std=c99")
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(CMAKE_EXE_LINKER_FLAGS_INIT "-Wl,--as-needed")
endif()
