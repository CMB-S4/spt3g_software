# - Find LAPACKE
# Find the native LAPACKE includes and library
#
# LAPACKE_INCLUDES - where to find fftw3.h
# LAPACKE_LIBRARIES - List of libraries when using LAPACKE.
# LAPACKE_FOUND - True if LAPACKE found.
if (LAPACKE_INCLUDES)
# Already in cache, be silent
set (LAPACKE_FIND_QUIETLY TRUE)
endif (LAPACKE_INCLUDES)
find_path (LAPACKE_INCLUDES lapacke.h)
find_library (LAPACKE_LIBRARIES NAMES lapacke)
# handle the QUIETLY and REQUIRED arguments and set LAPACKE_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (LAPACKE DEFAULT_MSG LAPACKE_LIBRARIES LAPACKE_INCLUDES)
mark_as_advanced (LAPACKE_LIBRARIES LAPACKE_INCLUDES)
