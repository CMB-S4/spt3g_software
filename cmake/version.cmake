set(CMAKE_SOURCE_DIR ${CMAKE_ARGV3})
set(CMAKE_BINARY_DIR ${CMAKE_ARGV4})

# Check VERSION file (sensible in exported source tree)
file(READ ${CMAKE_SOURCE_DIR}/VERSION GIT_VERSION)
string(REGEX REPLACE "\\$Version: (.*)\\$" "\\1" GIT_VERSION "${GIT_VERSION}")
string(REPLACE "$Version$" "" GIT_VERSION "${GIT_VERSION}")
string(STRIP "${GIT_VERSION}" GIT_VERSION)

# Cleanup on error
macro(cleanup)
    file(REMOVE ${CMAKE_BINARY_DIR}/cmake/Spt3gConfigVersion.cmake)
endmacro(cleanup)

if (NOT GIT_VERSION)
    find_package(Git QUIET)

    if (NOT GIT_EXECUTABLE)
        cleanup()
        return()
    endif()

    # Get version string from git
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_VERSION
        ERROR_VARIABLE GIT_VERSION_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )

    if (GIT_VERSION_ERROR)
        cleanup()
        return()
    endif()

endif()

# Extract numerical version, expect tag to be of the form v0.1.2
string(REGEX REPLACE "^v([0-9\\.]+).*" "\\1" VERSION "${GIT_VERSION}")
string(REGEX MATCH "^([0-9\\.]+)$" TEST_VERSION ${VERSION})

if (NOT TEST_VERSION)
    cleanup()
    return()
endif()

message(STATUS "Building version ${VERSION}")

# Populate the config file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/cmake/Spt3gConfigVersion.cmake"
    VERSION "${VERSION}"
    COMPATIBILITY AnyNewerVersion
    )
