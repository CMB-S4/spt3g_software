find_package(Git)

if (GIT_EXECUTABLE)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_VERSION
        ERROR_VARIABLE GIT_VERSION_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )

    if (GIT_VERSION_ERROR)
        message(WARNING "Cannot determine project version: ${GIT_VERSION_ERROR}")

    else()
        string(REGEX REPLACE "^v([0-9\\.]+).*" "\\1" VERSION "${GIT_VERSION}")

        if (NOT VERSION)
            message(WARNING "Incompatible git version string format: ${GIT_VERSION}")

        else()
            include(CMakePackageConfigHelpers)
            write_basic_package_version_file(
                "${CMAKE_BINARY_DIR}/cmake/Spt3gConfigVersion.cmake"
                VERSION "${VERSION}"
                COMPATIBILITY AnyNewerVersion
                )
        endif()

    endif()

endif()
