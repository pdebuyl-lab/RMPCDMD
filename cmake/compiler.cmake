set(CMAKE_BUILD_TYPE_INIT "Release")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_INIT "-Wall -Wextra -Wconversion -std=f2008 -pedantic")
endif()
