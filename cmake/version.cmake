FIND_PACKAGE(Git)

message("hop")
message(${CMAKE_SOURCE_DIR}/..)

IF(GIT_EXECUTABLE)
  execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --dirty OUTPUT_VARIABLE RMPCDMD_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE)
  configure_file(${CMAKE_SOURCE_DIR}/../src/rmpcdmd_module.f90.in ${CMAKE_BINARY_DIR}/rmpcdmd_module.f90 @ONLY)
else()
  message(FATAL_ERROR "git not found")
endif(GIT_EXECUTABLE)
