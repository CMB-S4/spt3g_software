# Locate Python

if (BUILD_WHEEL)
	find_package(Python COMPONENTS Interpreter REQUIRED)
elseif(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.12)
	find_package(Python COMPONENTS Interpreter Development)
else()
	find_package(PythonInterp)
	find_package(PythonLibs ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
	string(REGEX REPLACE ".*libpython([0-9])\\.[0-9]+.*\\..*" "\\1" Python_VERSION_MAJOR ${PYTHON_LIBRARIES})
	string(REGEX REPLACE ".*libpython[0-9]\\.([0-9]+).*\\..*" "\\1" Python_VERSION_MINOR ${PYTHON_LIBRARIES})
	set(Python_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH})
	set(Python_LIBRARIES ${PYTHON_LIBRARIES})
	set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
endif()

# look for numpy
execute_process(COMMAND ${Python_EXECUTABLE} -c "import numpy"
	RESULT_VARIABLE NUMPY_FOUND ERROR_QUIET)
if(NUMPY_FOUND EQUAL 0)
	set(Python_NumPy_FOUND TRUE CACHE BOOL "Numpy found successfully" FORCE)
else(NUMPY_FOUND EQUAL 0)
	set(Python_NumPy_FOUND FALSE CACHE BOOL "Numpy found successfully" FORCE)
endif(NUMPY_FOUND EQUAL 0)

if(Python_NumPy_FOUND)
	execute_process(COMMAND ${Python_EXECUTABLE} -c
		"import numpy; print(numpy.get_include())"
		OUTPUT_VARIABLE _NUMPY_INCLUDE_DIR
		OUTPUT_STRIP_TRAILING_WHITESPACE)

	find_path(NUMPY_INCLUDE_DIR NAMES numpy/ndarrayobject.h HINTS ${_NUMPY_INCLUDE_DIR})
	set(Python_NumPy_INCLUDE_DIRS ${NUMPY_INCLUDE_DIR} CACHE STRING "Numpy inc directory")
	message(STATUS "Found NumPy: ${Python_NumPy_INCLUDE_DIRS}")
else(Python_NumPy_FOUND)
	message(STATUS "Found NumPy: NOT FOUND")
endif(Python_NumPy_FOUND)

## look for scipy
execute_process(COMMAND ${Python_EXECUTABLE} -c "import scipy"
	RESULT_VARIABLE SCIPY_FOUND OUTPUT_QUIET ERROR_QUIET)
if(SCIPY_FOUND EQUAL 0)
	set(Python_SciPy_FOUND TRUE CACHE BOOL "SciPy Found")
	message(STATUS "Found SciPy: FOUND")
else()
	set(Python_SciPy_FOUND FALSE CACHE BOOL "SciPy Found")
	message(STATUS "Found SciPy: NOT FOUND")
endif()

# suppress configuration warnings in newer cmake / boost versions
if(NOT DEFINED CMAKE_FIND_PACKAGE_PREFER_CONFIG)
	set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)
endif()

if(NOT DEFINED Boost_PYTHON_TYPE)
	set(Boost_PYTHON_TYPE python)

	# Hack for some old Boost CMake modules
	set(_Boost_PYTHON${Python_VERSION_MAJOR}_HEADERS "boost/python.hpp")
	set(_Boost_PYTHON${Python_VERSION_MAJOR}${Python_VERSION_MINOR}_HEADERS "boost/python.hpp")

	find_package(Boost QUIET COMPONENTS python${Python_VERSION_MAJOR}${Python_VERSION_MINOR})
	if (${Boost_PYTHON${Python_VERSION_MAJOR}${Python_VERSION_MINOR}_FOUND})
		set(Boost_PYTHON_TYPE python${Python_VERSION_MAJOR}${Python_VERSION_MINOR})
	else()
		find_package(Boost QUIET COMPONENTS python${Python_VERSION_MAJOR})
		if (${Boost_PYTHON${Python_VERSION_MAJOR}_FOUND})
			set(Boost_PYTHON_TYPE python${Python_VERSION_MAJOR})
		endif()
	endif()
endif()

find_package(Boost COMPONENTS system iostreams filesystem ${Boost_PYTHON_TYPE} REQUIRED)
message(STATUS "Found Boost: ${Boost_INCLUDE_DIR} (found version \"${Boost_VERSION}\")")
