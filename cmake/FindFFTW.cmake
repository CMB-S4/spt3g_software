# - Find FFTW
# Find the native FFTW includes and library
#
# FFTW_INCLUDE_DIRS - where to find fftw3.h
# FFTW_LIBRARIES - List of libraries when using FFTW.
# FFTW_FOUND - True if FFTW found.
if (FFTW_INCLUDE_DIRS)
# Already in cache, be silent
set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDE_DIRS)
find_path (FFTW_INCLUDE_DIRS fftw3.h HINTS ENV FFTW_INC)
find_library (FFTW_LIBRARIES NAMES fftw3 HINTS ENV FFTW_DIR)


#handles finding thread implementations
find_library (FFTW_THREADS_LIBRARY NAMES fftw3_threads HINTS ENV FFTW_DIR)

include (FindPackageHandleStandardArgs)
set(FPHSA_NAME_MISMATCHED TRUE)
find_package_handle_standard_args(FFTW_THREADS DEFAULT_MSG FFTW_THREADS_LIBRARY)
mark_as_advanced(FFTW_THREADS_LIBRARY)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDE_DIRS)
mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDE_DIRS)
