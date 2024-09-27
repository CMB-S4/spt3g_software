cmake_minimum_required(VERSION 3.12)
# Don't warn about removal of FindBoost in cmake 3.30+
if(POLICY CMP0167)
	cmake_policy(SET CMP0167 NEW)
endif()

project(spt3g)
enable_testing()

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# Use cmake directory for packages
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{SITE_CMAKE_DIR} "${CMAKE_SOURCE_DIR}/cmake")

# Work around stupid broken Red Hat systems
set(CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "")

# Require C++ 14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Convenience variables
set(SPT3G_MODULE_DIR ${CMAKE_BINARY_DIR}/spt3g)
set(SPT3G_INCLUDE_INSTALL_DIR "include/spt3g")

include(Spt3gIncludes)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Raise every warning by default
# (use target-specific options to disable particular warnings)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find all the Boost and Python libraries
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

if(DEFINED ENV{CIBUILDWHEEL})
	# Can't compile binaries against libpython in CI
	find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
else()
	find_package(Python COMPONENTS Interpreter Development REQUIRED)
endif()

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(Boost_PYTHON_VERSION ${Python_VERSION})
find_package(Boost COMPONENTS system iostreams filesystem python REQUIRED)

# Interface library for flags and library dependencies
add_library(spt3g INTERFACE)

# get version info
include(Spt3gVersion)
if (SPT3G_VERSION)
	message(STATUS "Building SPT3G software version ${SPT3G_VERSION}")
	target_compile_definitions(spt3g INTERFACE -DSPT3G_VERSION="${SPT3G_VERSION}")
endif()

target_compile_features(spt3g INTERFACE cxx_std_11)
target_include_directories(spt3g INTERFACE ${Boost_INCLUDE_DIR})
target_link_libraries(spt3g INTERFACE Threads::Threads ${Boost_LIBRARIES})
target_include_directories(spt3g INTERFACE $<INSTALL_INTERFACE:${SPT3G_INCLUDE_INSTALL_DIR}>)

# Python bindings
if(TARGET Python::Module)
	target_link_libraries(spt3g INTERFACE Python::Module)
else()
	target_include_directories(spt3g INTERFACE ${Python_INCLUDE_DIRS})
	if(APPLE)
		target_link_options(spt3g INTERFACE "LINKER:-undefined,dynamic_lookup")
	endif()
endif()

# Work around yet more bugs in GCC 4.4, this time with C++ 11 support
# Also increase maximum number of arguments in python bindings
target_compile_definitions(spt3g INTERFACE -DBOOST_PYTHON_MAX_ARITY=20
	-DBOOST_ALLOW_DEPRECATED_HEADERS)
if(Boost_VERSION VERSION_LESS 1.77)
	target_compile_definitions(spt3g INTERFACE -DBOOST_BIND_GLOBAL_PLACEHOLDERS)
endif()
if(Boost_VERSION VERSION_LESS 1.82)
	target_compile_definitions(spt3g INTERFACE -DBOOST_NO_CXX11_SMART_PTR=1)
endif()

# Shell script to set environment variables
configure_file(${CMAKE_SOURCE_DIR}/cmake/env-shell.sh.in ${CMAKE_BINARY_DIR}/env-shell.sh @ONLY)

# Set up python importability
execute_process(COMMAND mkdir -p ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
execute_process(COMMAND mkdir -p ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
execute_process(COMMAND mkdir -p ${SPT3G_MODULE_DIR})
execute_process(COMMAND ln -fsn ${CMAKE_SOURCE_DIR}/cmake/package/__init__.py ${SPT3G_MODULE_DIR}/__init__.py)

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

include(GNUInstallDirs)
if(APPLE)
       # See: https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling
	set(CMAKE_MACOSX_RPATH TRUE)
	set(CMAKE_SKIP_BUILD_RPATH FALSE)
	set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
	set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
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
if(NOT PIP_SPT3G_VERSION_FILE)
	add_custom_target(version ALL
		COMMAND sh ${CMAKE_SOURCE_DIR}/cmake/getvers.sh ${CMAKE_SOURCE_DIR} ${SPT3G_MODULE_DIR}/version.py
		BYPRODUCTS ${SPT3G_MODULE_DIR}/version.py
		COMMENT "Regenerating VCS version information"
	)
endif()

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

set(INCLUDE_INSTALL_DIR ${SPT3G_INCLUDE_INSTALL_DIR} )
set(LIB_INSTALL_DIR lib/ )
set(CMAKE_INSTALL_DIR lib/cmake/${PROJECT_NAME} )
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
	"${CMAKE_SOURCE_DIR}/cmake/FindFLAC.cmake"
	"${CMAKE_SOURCE_DIR}/cmake/FindNetCDF.cmake"
	"${CMAKE_SOURCE_DIR}/cmake/FindOgg.cmake"
	DESTINATION ${CMAKE_INSTALL_DIR})
if(SPT3G_VERSION)
	install(FILES "${CMAKE_BINARY_DIR}/cmake/Spt3gConfigVersion.cmake" DESTINATION ${CMAKE_INSTALL_DIR})
endif()

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
