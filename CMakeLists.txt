cmake_minimum_required(VERSION 2.6)
if(POLICY CMP0060) # Avoid cmake misfeature that can link the wrong libraries
	cmake_policy(SET CMP0060 NEW)
endif()
cmake_policy(SET CMP0012 NEW) # Allow use of true in boolean expressions

project(spt3g)
enable_testing()

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spt3g)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

# Work around stupid broken Red Hat systems
set(CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
	message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler (GCC > 4.7, Clang).")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unknown-warning-option -Wno-unused -Wno-unused-result -Wno-deprecated-register -Wno-sign-compare -Werror -Wno-constant-conversion -Wno-self-assign-overloaded")

# Fix bugs in GCC 4.4's strict aliasing code by turning it off
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing")

# Work around yet more bugs in GCC 4.4, this time with C++ 11 support
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_NO_CXX11_SMART_PTR=1")

# Boost bits we need
find_package(PythonInterp)
find_package(PythonLibs ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(Boost_NO_BOOST_CMAKE ON)

if(NOT DEFINED Boost_PYTHON_TYPE)
	set(Boost_PYTHON_TYPE python)

	# Detect attempts to be clever with the naming of the Boost Python library
	string(REGEX REPLACE ".*libpython([0-9])\\.[0-9]+.*\\..*" "\\1" PYTHONMAJORVER ${PYTHON_LIBRARIES})
	string(REGEX REPLACE ".*libpython[0-9]\\.([0-9]+).*\\..*" "\\1" PYTHONMINORVER ${PYTHON_LIBRARIES})

	# Hack for some old Boost CMake modules
	set(_Boost_PYTHON${PYTHONMAJORVER}_HEADERS "boost/python.hpp")
	set(_Boost_PYTHON${PYTHONMAJORVER}${PYTHONMINORVER}_HEADERS "boost/python.hpp")

	find_package(Boost QUIET COMPONENTS python${PYTHONMAJORVER}${PYTHONMINORVER})
	if (${Boost_PYTHON${PYTHONMAJORVER}${PYTHONMINORVER}_FOUND})
		set(Boost_PYTHON_TYPE python${PYTHONMAJORVER}${PYTHONMINORVER})
	else()
		find_package(Boost QUIET COMPONENTS python${PYTHONMAJORVER})
		if (${Boost_PYTHON${PYTHONMAJORVER}_FOUND})
			set(Boost_PYTHON_TYPE python${PYTHONMAJORVER})
		endif()
	endif()
endif()

find_package(Boost COMPONENTS system iostreams filesystem ${Boost_PYTHON_TYPE} REQUIRED)
include_directories(${Boost_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH})
if(Boost_VERSION EQUAL 104700 OR Boost_VERSION GREATER 104700)
	set(Boost_GEOMETRY_AVAIL ON)
	message(STATUS "Boost Geometry Available")
endif()


# Make python import more idiomatic
set(CMAKE_SHARED_LIBRARY_PREFIX "")

# Shell script to set environment variables
configure_file(${CMAKE_SOURCE_DIR}/cmake/env-shell.sh.in ${CMAKE_BINARY_DIR}/env-shell.sh @ONLY)

# Set up python importability
execute_process(COMMAND mkdir -p ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
execute_process(COMMAND mkdir -p ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
execute_process(COMMAND ln -fsn ${CMAKE_SOURCE_DIR}/cmake/init.py ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/__init__.py)

# Use cmake directory for packages
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{SITE_CMAKE_DIR} "${CMAKE_SOURCE_DIR}/cmake")

# Frequently used external code
set(OpenMP_FIND_QUIETLY TRUE)
find_package(GSL)
find_package(OpenMP)
find_package(FFTW)
find_package(HDF5)
find_package(MPI)

# Convenience macros
macro(link_python_dir)
	execute_process(COMMAND mkdir -p ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
	if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/python)
		execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR}/python ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${PROJECT})
	else(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/python)
		execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${PROJECT})
	endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/python)
endmacro(link_python_dir)
macro(add_spt3g_test test_name)
	add_test(${PROJECT}/${test_name} ${CMAKE_CURRENT_SOURCE_DIR}/tests/${test_name}.py)

    set (extra_macro_args ${ARGN})
    list(LENGTH extra_macro_args num_extra_args)
    if (${num_extra_args} GREATER 0)
        list(GET extra_macro_args 0 test_labels)
		set_tests_properties(${PROJECT}/${test_name} PROPERTIES LABELS ${test_labels})
    endif ()
endmacro(add_spt3g_test test_name)

# Find all sub-projects
file(GLOB cmake_projects RELATIVE ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/*/CMakeLists.txt)
foreach(d ${cmake_projects})
	get_filename_component(proj ${d} PATH)
	set(SUBDIRS ${SUBDIRS} ${proj})
	get_filename_component(pname ${proj} NAME_WE)
	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${pname}/include)
endforeach(d ${cmake_projects})

# Make core includes available without directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/core/include/core)

list(SORT SUBDIRS)
foreach(subdir ${SUBDIRS})
	get_filename_component(pname ${subdir} NAME_WE)
	message(STATUS "+ ${pname}")
	set(PROJECT ${pname})
	add_subdirectory(${CMAKE_SOURCE_DIR}/${pname})
endforeach(subdir ${SUBDIRS})

# Add fetching of test data
add_custom_target(testdata COMMAND rsync -vrlpt --delete rsync://bolo.berkeley.edu/testdata ${CMAKE_BINARY_DIR}/testdata COMMENT "Rsyncing test data from bolo.berkeley.edu")

add_custom_target(fasttest COMMAND ctest -LE SLOWTEST COMMENT "Running fast test suite")

# Custom things related to testing
configure_file(${CMAKE_SOURCE_DIR}/cmake/test_env.py.in ${CMAKE_BINARY_DIR}/test_env.py @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/cmake/CTestCustom.cmake.in ${CMAKE_BINARY_DIR}/CTestCustom.cmake @ONLY)

# Target for version string
add_custom_target(version ALL
	COMMAND sh ${CMAKE_SOURCE_DIR}/cmake/getvers.sh ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR}/spt3g/__version__.py
	BYPRODUCTS ${CMAKE_BINARY_DIR}/spt3g/__version__.py
	COMMENT "Regenerating VCS version information"
)

# Add mechanism to make a tarball for the grid
add_custom_target(tarball
	COMMAND rm -rf ${CMAKE_BINARY_DIR}/tarball ${CMAKE_BINARY_DIR}/spt3g.tgz
	COMMAND mkdir ${CMAKE_BINARY_DIR}/tarball
	COMMAND ${CMAKE_SOURCE_DIR}/cmake/copy_build_dir.sh ${CMAKE_BINARY_DIR}
	COMMAND ${PYTHON_EXECUTABLE} -m compileall -fq ${CMAKE_BINARY_DIR}/tarball
	COMMAND tar -C ${CMAKE_BINARY_DIR}/tarball/ -czf ${CMAKE_BINARY_DIR}/spt3g.tgz .
	COMMENT "Build a tarball to run on the grid in spt3g.tgz")

# Add target to generate documentation
add_custom_target(docs
	COMMAND cmake -P ${CMAKE_SOURCE_DIR}/cmake/gen_per_module_docs.cmake ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR}
	COMMAND mkdir -p ${CMAKE_BINARY_DIR}/docs
	COMMAND which sphinx-build && sphinx-build -b html ${CMAKE_SOURCE_DIR}/doc ${CMAKE_BINARY_DIR}/docs || echo "Error: sphinx-build not found. Could not generate HTML docs."
	COMMENT "Generate HTML documentation")
