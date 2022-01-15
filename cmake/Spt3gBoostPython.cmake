# Locate Python
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.14)
	find_package(Python COMPONENTS Interpreter Development NumPy)
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
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.14)
	set(NUMPY_FOUND ${Python_NumPy_FOUND} CACHE BOOL "Numpy found successfully")
	set(NUMPY_INCLUDE_DIR ${Python_NumPy_INCLUDE_DIRS})
	message(STATUS "Found NumPy: ${Python_NumPy_INCLUDE_DIRS}")
else(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.14)
	# Old, crummy cmake -- try our best. Use better cmake for special cases.

	execute_process(COMMAND ${Python_EXECUTABLE} -c "import numpy"
		RESULT_VARIABLE NUMPY_FOUND)
	if(NUMPY_FOUND EQUAL 0)
		set(Python_NumPy_FOUND TRUE CACHE BOOL "Numpy found successfully")
	else(NUMPY_FOUND EQUAL 0)
		set(Python_NumPy_FOUND FALSE CACHE BOOL "Numpy found successfully")
	endif(NUMPY_FOUND EQUAL 0)

	execute_process(COMMAND ${Python_EXECUTABLE} -c
		"import numpy; print(numpy.get_include())"
		OUTPUT_VARIABLE _NUMPY_INCLUDE_DIR
		OUTPUT_STRIP_TRAILING_WHITESPACE)

	find_path(NUMPY_INCLUDE_DIR NAMES numpy/ndarrayobject.h HINTS ${_NUMPY_INCLUDE_DIR})
	set(Python_NumPy_INCLUDE_DIRS ${NUMPY_INCLUDE_DIR} CACHE STRING "Numpy inc directory")
	message(STATUS "Found NumPy: ${Python_NumPy_INCLUDE_DIRS}")
endif()

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
set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)

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
