set(CMAKE_SOURCE_DIR ${CMAKE_ARGV3})
set(CMAKE_BINARY_DIR ${CMAKE_ARGV4})

file(READ ${CMAKE_SOURCE_DIR}/VERSION GIT_VERSION)
string(STRIP "${GIT_VERSION}" GIT_VERSION)

if (NOT GIT_VERSION)
    find_package(Git QUIET)

    if (NOT GIT_EXECUTABLE)
        message(WARNING "Cannot determine project version: Missing git executable")
        return()
    endif()

    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_VERSION
        ERROR_VARIABLE GIT_VERSION_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )

    if (GIT_VERSION_ERROR)
        message(WARNING "Cannot determine project version: ${GIT_VERSION_ERROR}")
        return()
    endif()

endif()

string(REGEX REPLACE "^v([0-9\\.]+).*" "\\1" VERSION "${GIT_VERSION}")
string(REGEX MATCH "[0-9\\.]+" TEST_VERSION ${VERSION})

if (NOT TEST_VERSION)
    message(WARNING "Incompatible git version string: ${GIT_VERSION}")
    return()
endif()

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/cmake/Spt3gConfigVersion.cmake"
    VERSION "${VERSION}"
    COMPATIBILITY AnyNewerVersion
    )
