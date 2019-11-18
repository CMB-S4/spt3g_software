if (NOT CMAKE_PARENT_LIST_FILE)
	set(CMAKE_SOURCE_DIR ${CMAKE_ARGV3})
	set(CMAKE_BINARY_DIR ${CMAKE_ARGV4})
endif()

function(get_spt3g_version)
	# Check VERSION file (sensible in exported source tree)
	file(READ ${CMAKE_SOURCE_DIR}/VERSION GIT_VERSION)
	string(REGEX REPLACE "\\$Version: (.*)\\$" "\\1" GIT_VERSION "${GIT_VERSION}")
	string(REPLACE "$Version$" "" GIT_VERSION "${GIT_VERSION}")
	string(STRIP "${GIT_VERSION}" GIT_VERSION)

	if (NOT GIT_VERSION)
		find_package(Git QUIET)

		if (NOT GIT_EXECUTABLE)
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
			return()
		endif()

	endif()

	# Extract numerical version, expect tag to be of the form v0.1.2
	string(REGEX REPLACE "^v([0-9\\.]+).*" "\\1" VERSION "${GIT_VERSION}")
	string(REGEX MATCH "^([0-9\\.]+)$" TEST_VERSION ${VERSION})

	if (NOT TEST_VERSION)
		return()
	endif()

	set(SPT3G_VERSION ${VERSION} PARENT_SCOPE)
endfunction(get_spt3g_version)

if (NOT CMAKE_PARENT_LIST_FILE)
	get_spt3g_version()
	if (NOT SPT3G_VERSION)
		file(REMOVE ${CMAKE_BINARY_DIR}/cmake/Spt3gConfigVersion.cmake)
		return()
	endif()

	message(STATUS "Building version ${SPT3G_VERSION}")

	# Populate the config file
	include(CMakePackageConfigHelpers)
	write_basic_package_version_file(
		"${CMAKE_BINARY_DIR}/cmake/Spt3gConfigVersion.cmake"
		VERSION "${SPT3G_VERSION}"
		COMPATIBILITY AnyNewerVersion
		)
endif()
