IF(NOT HDF5_ROOT)
MESSAGE(FATAL_ERROR "HDF5_ROOT undefined, stopping")
ENDIF(NOT HDF5_ROOT)
IF(NOT Z_ROOT)
SET(Z_ROOT ${HDF5_ROOT})
ENDIF(NOT Z_ROOT)
IF(NOT SZ_ROOT)
SET(SZ_ROOT ${HDF5_ROOT})
ENDIF(NOT SZ_ROOT)

FIND_LIBRARY(HDF5LIB hdf5 HINTS ${HDF5_ROOT} PATH_SUFFIXES lib include NO_DEFAULT_PATH)
FIND_LIBRARY(HDF5FORTRANLIB hdf5_fortran HINTS ${HDF5_ROOT} PATH_SUFFIXES lib include NO_DEFAULT_PATH)
FIND_LIBRARY(ZLIB z HINTS ${Z_ROOT} PATH_SUFFIXES lib include NO_DEFAULT_PATH)
FIND_LIBRARY(SZLIB sz HINTS ${SZ_ROOT} PATH_SUFFIXES lib include NO_DEFAULT_PATH)
FIND_PATH(HDF5_INC hdf5.mod HINTS ${HDF5_ROOT} PATH_SUFFIXES lib include NO_DEFAULT_PATH)
