cmake_minimum_required(VERSION 3.5)
if(POLICY CMP0060) # Avoid cmake misfeature that can link the wrong libraries
	cmake_policy(SET CMP0060 NEW)
endif()
cmake_policy(SET CMP0012 NEW) # Allow use of true in boolean expressions
if(POLICY CMP0042) # Enable RPATH on OSX
	cmake_policy(SET CMP0042 NEW)
endif()
if(POLICY CMP0057) # Allow IN_LIST operator in if() statements
	cmake_policy(SET CMP0057 NEW)
endif()

project(spt3g)
enable_testing()

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# Work around stupid broken Red Hat systems
set(CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "")

# Require C++ 14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler and include configurations
set(SPT3G_SOURCE_DIR ${CMAKE_SOURCE_DIR})
set(SPT3G_BUILD_DIR ${CMAKE_BINARY_DIR})

set(SPT3G_INCLUDE_INSTALL_DIR "include/spt3g")

configure_file(${CMAKE_SOURCE_DIR}/cmake/Spt3gIncludes.cmake.in ${CMAKE_BINARY_DIR}/cmake/Spt3gIncludes.cmake @ONLY)
include(${CMAKE_BINARY_DIR}/cmake/Spt3gIncludes.cmake)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${SPT3G_LIBRARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SPT3G_RUNTIME_DIR})

# Sneak in a namespace-like prefix for all libraries
set(CMAKE_SHARED_LIBRARY_PREFIX "libspt3g-")

# Raise errors on every warning by default
# (use target-specific options to disable particular warnings)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")

# Interface library for flags and library dependencies
add_library(spt3g INTERFACE)

# get version info
include(Spt3gVersion)
get_spt3g_version()
if (SPT3G_VERSION)
	message(STATUS "Building SPT3G software version ${SPT3G_VERSION}")
	target_compile_definitions(spt3g INTERFACE -DSPT3G_VERSION="${SPT3G_VERSION}")
endif()

# Ignore known warnings
target_compile_options(spt3g INTERFACE -Wno-unknown-warning-option -Wno-unused -Wno-unused-result -Wno-deprecated-register -Wno-sign-compare -Wno-constant-conversion -Wno-self-assign-overloaded -Wno-deprecated-declarations -Wno-maybe-uninitialized)

# Fix bugs in GCC 4.4's strict aliasing code by turning it off
target_compile_options(spt3g INTERFACE -fno-strict-aliasing)

# Boost bits we need
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(Boost_NO_BOOST_CMAKE ON)

# Find all the Boost and Python libraries
include(Spt3gBoostPython)

# Upgrade C++ standard to C++14 for Boost >= 1.80, which requires it
if(Boost_VERSION VERSION_GREATER 1.79)
	set(CMAKE_CXX_STANDARD 14)
endif()

target_include_directories(spt3g INTERFACE ${Boost_INCLUDE_DIR} ${Python_INCLUDE_DIRS})
target_link_libraries(spt3g INTERFACE pthread ${Boost_LIBRARIES} ${Python_LIBRARIES})

# Work around yet more bugs in GCC 4.4, this time with C++ 11 support
# Also increase maximum number of arguments in python bindings
target_compile_definitions(spt3g INTERFACE -DBOOST_NO_CXX11_SMART_PTR=1
    -DBOOST_PYTHON_MAX_ARITY=20 -DBOOST_BIND_GLOBAL_PLACEHOLDERS
    -DBOOST_ALLOW_DEPRECATED_HEADERS)

# Shell script to set environment variables
configure_file(${CMAKE_SOURCE_DIR}/cmake/env-shell.sh.in ${CMAKE_BINARY_DIR}/env-shell.sh @ONLY)

# Set up python importability
execute_process(COMMAND mkdir -p ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
execute_process(COMMAND mkdir -p ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
execute_process(COMMAND ln -fsn ${CMAKE_SOURCE_DIR}/cmake/init.py ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/__init__.py)

set(BUILD_PROJECTS "${BUILD_PROJECTS}" CACHE STRING "The subset of available projects to actually build")
if(NOT "${BUILD_PROJECTS}" STREQUAL "")
	set(USE_PROJECT_LIST "TRUE" CACHE BOOL "Whether to build only the subset of projects sepcified by BUILD_PROJECTS")
else(NOT "${BUILD_PROJECTS}" STREQUAL "")
	set(USE_PROJECT_LIST "FALSE" CACHE BOOL "Whether to build only the subset of projects sepcified by BUILD_PROJECTS")
endif(NOT "${BUILD_PROJECTS}" STREQUAL "")

list(FIND BUILD_PROJECTS core WILL_BUILD_CORE)
if(${USE_PROJECT_LIST} AND WILL_BUILD_CORE EQUAL -1)
	message(STATUS "Automatically adding the core project to the manually specified BUILD_PROJECTS list")
	list(APPEND BUILD_PROJECTS core)
endif(${USE_PROJECT_LIST} AND WILL_BUILD_CORE EQUAL -1)

# If not set by the user, ask python where it keeps modules
if(NOT PYTHON_MODULE_DIR)
	EXECUTE_PROCESS (COMMAND ${Python_EXECUTABLE} -c
		"import sysconfig; print(sysconfig.get_path(\"platlib\"))"
		OUTPUT_VARIABLE PYTHON_MODULE_DIR
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
endif(NOT PYTHON_MODULE_DIR)
message(STATUS "Python modules will be installed to ${PYTHON_MODULE_DIR}")

include(GNUInstallDirs)
if(APPLE)
       # See: https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling
	set(CMAKE_MACOSX_RPATH TRUE)
	set(CMAKE_SKIP_BUILD_RPATH FALSE)
	set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
	set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${PYTHON_MODULE_DIR}/spt3g")
	set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif(APPLE)

# Find all sub-projects
file(GLOB cmake_projects RELATIVE ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/*/CMakeLists.txt)
foreach(d ${cmake_projects})
	get_filename_component(proj ${d} PATH)
	set(SUBDIRS ${SUBDIRS} ${proj})
endforeach(d ${cmake_projects})

set(SPT3G_LIBRARIES spt3g)
set(SPT3G_PROGRAMS "")
set(SPT3G_PYTHON_DIRS "")
list(SORT SUBDIRS)
foreach(subdir ${SUBDIRS})
	get_filename_component(pname ${subdir} NAME_WE)

	list(FIND BUILD_PROJECTS "${pname}" SHOULD_BUILD_PROJECT)
	if(${USE_PROJECT_LIST} AND SHOULD_BUILD_PROJECT EQUAL -1)
		continue() # if we're using the project list but this one isn't in it, skip it
	endif(${USE_PROJECT_LIST} AND SHOULD_BUILD_PROJECT EQUAL -1)

	message(STATUS "+ ${pname}")
	set(PROJECT ${pname})
	add_subdirectory(${CMAKE_SOURCE_DIR}/${pname})
endforeach(subdir ${SUBDIRS})

# export configuration files for use in other projects
export(TARGETS ${SPT3G_LIBRARIES} NAMESPACE spt3g:: FILE ${CMAKE_BINARY_DIR}/cmake/Spt3gTargets.cmake)
configure_file(${CMAKE_SOURCE_DIR}/cmake/Spt3gConfig.cmake.in ${CMAKE_BINARY_DIR}/cmake/Spt3gConfig.cmake @ONLY)

# Custom things related to testing
add_custom_target(fasttest COMMAND ctest -LE SLOWTEST COMMENT "Running fast test suite")
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.17)
	list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure")
endif(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.17)

# Target for version string
add_custom_target(version ALL
	COMMAND sh ${CMAKE_SOURCE_DIR}/cmake/getvers.sh ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR}
	COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/cmake/Spt3gVersion.cmake ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR}
	BYPRODUCTS ${CMAKE_BINARY_DIR}/spt3g/version.py ${CMAKE_BINARY_DIR}/pyproject.toml ${CMAKE_BINARY_DIR}/cmake/Spt3gConfigVersion.cmake
	COMMENT "Regenerating VCS version information"
)

# Add mechanism to make a tarball for the grid
add_custom_target(tarball
	COMMAND rm -rf ${CMAKE_BINARY_DIR}/tarball ${CMAKE_BINARY_DIR}/spt3g.tgz
	COMMAND mkdir ${CMAKE_BINARY_DIR}/tarball
	COMMAND ${CMAKE_SOURCE_DIR}/cmake/copy_build_dir.sh ${CMAKE_BINARY_DIR}
	COMMAND ${Python_EXECUTABLE} -m compileall -fq ${CMAKE_BINARY_DIR}/tarball
	COMMAND tar -C ${CMAKE_BINARY_DIR}/tarball/ -czf ${CMAKE_BINARY_DIR}/spt3g.tgz .
	COMMENT "Build a tarball to run on the grid in spt3g.tgz")

# Add target to generate documentation
add_custom_target(docs
	COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/cmake/gen_per_module_docs.cmake ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} ${Python_EXECUTABLE}
	COMMAND mkdir -p ${CMAKE_BINARY_DIR}/docs
	COMMAND which sphinx-build && sphinx-build -b html ${CMAKE_SOURCE_DIR}/doc ${CMAKE_BINARY_DIR}/docs || echo "Error: sphinx-build not found. Could not generate HTML docs."
	COMMENT "Generate HTML documentation")

# Set up installation

# Export the phony library target
INSTALL(TARGETS spt3g EXPORT ${PROJECT_NAME}Config)

# The python module
# First, one-off files
INSTALL(FILES ${CMAKE_SOURCE_DIR}/cmake/init.py DESTINATION ${PYTHON_MODULE_DIR}/spt3g RENAME __init__.py)
INSTALL(FILES ${CMAKE_BINARY_DIR}/spt3g/version.py DESTINATION ${PYTHON_MODULE_DIR}/spt3g)

# Python scripts
foreach(dir ${SPT3G_PYTHON_DIRS})
	get_filename_component(dir_name "${dir}" NAME)
	if("${dir_name}" STREQUAL "python")
		# If the directory name is "python", need to walk up the path to find out what it should be called.
		get_filename_component(dir_name "${dir}" DIRECTORY)
		get_filename_component(lib "${dir_name}" NAME)
		# The INSTALL conmmand does not support the RENAME option in directory mode, which makes this a pain.
		file(GLOB_RECURSE project_py_files LIST_DIRECTORIES false RELATIVE ${dir} ${dir}/*)
		foreach(file ${project_py_files})
			get_filename_component(subdir "${file}" DIRECTORY)
			if(file MATCHES ".*pyc" OR file MATCHES ".*pyo")
				continue() # skip existing compiled python files
			endif(file MATCHES ".*pyc" OR file MATCHES ".*pyo")
			install(FILES "${dir}/${file}" DESTINATION ${PYTHON_MODULE_DIR}/spt3g/${lib}/${subdir})
		endforeach(file ${project_py_files})
	else("${dir_name}" STREQUAL "python")
		install(DIRECTORY ${dir} DESTINATION ${PYTHON_MODULE_DIR}/spt3g/
			PATTERN "*.pyc" EXCLUDE # skip existing compiled python files
			PATTERN "*.pyo" EXCLUDE)
	endif("${dir_name}" STREQUAL "python")
endforeach(dir SPT3G_PYTHON_DIRS)
# Ensure that python scripts are precompiled, in case they will be executed from a read-only perspective
install(CODE "MESSAGE(STATUS \"Pre-compiling python scripts in ${PYTHON_MODULE_DIR}/spt3g/\")
EXECUTE_PROCESS(COMMAND ${Python_EXECUTABLE} -O -m compileall -fq ${PYTHON_MODULE_DIR}/spt3g/)")

# The exectutables
foreach(program ${SPT3G_PROGRAMS})
	# programs are likely to by symlinks, so resolve those before installing
	if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.19)
		file(REAL_PATH "${CMAKE_BINARY_DIR}/bin/${program}" true_program)
	else(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.19)
		get_filename_component(true_program "${CMAKE_BINARY_DIR}/bin/${program}" REALPATH)
	endif(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.19)
	install(PROGRAMS "${true_program}"
		DESTINATION "${CMAKE_INSTALL_BINDIR}"
		RENAME "${program}"
	)
endforeach(program ${SPT3G_PROGRAMS})

# The CMake package definitions
include(CMakePackageConfigHelpers)

set(INCLUDE_INSTALL_DIR include/spt3g )
set(LIB_INSTALL_DIR lib/ )
set(CMAKE_INSTALL_DIR share/${PROJECT_NAME}/cmake )
configure_package_config_file(${CMAKE_SOURCE_DIR}/cmake/Spt3gConfig.cmake.in
	${CMAKE_BINARY_DIR}/cmake-install/Spt3gConfig.cmake
	INSTALL_DESTINATION ${CMAKE_INSTALL_DIR}
	PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR CMAKE_INSTALL_DIR
)

install(EXPORT ${PROJECT_NAME}Config
	DESTINATION ${CMAKE_INSTALL_DIR}
	NAMESPACE spt3g::
	FILE Spt3gTargets.cmake)
install(FILES
	"${CMAKE_BINARY_DIR}/cmake-install/Spt3gConfig.cmake"
	"${CMAKE_BINARY_DIR}/cmake/Spt3gIncludes.cmake"
	DESTINATION ${CMAKE_INSTALL_DIR})

# The header files
foreach(lib ${SPT3G_LIBRARIES})
	if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${lib}/include/)
		file(GLOB INCLUDE_ITEMS LIST_DIRECTORIES TRUE
			RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/${lib}/include/
			CONFIGURE_DEPENDS
			"${CMAKE_CURRENT_SOURCE_DIR}/${lib}/include/*")
		foreach(item ${INCLUDE_ITEMS})
			install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${lib}/include/${item}
				DESTINATION ${INCLUDE_INSTALL_DIR}
			)
		endforeach(item ${INCLUDE_ITEMS})
	endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${lib}/include/)
endforeach(lib ${SPT3G_LIBRARIES})

# cmake doesn't generally have headers, but we want to provide the test header
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/G3Test.h DESTINATION ${INCLUDE_INSTALL_DIR})

# Target for uninstalling
configure_file(
	"${CMAKE_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
	"${CMAKE_BINARY_DIR}/cmake/cmake_uninstall_cmake"
	IMMEDIATE @ONLY)
add_custom_target(uninstall "${CMAKE_COMMAND}" -P "${CMAKE_BINARY_DIR}/cmake/cmake_uninstall_cmake")

if(NOT Python_FOUND)
	message(FATAL_ERROR "Unable to find Python; build cannot proceed")
endif(NOT Python_FOUND)
