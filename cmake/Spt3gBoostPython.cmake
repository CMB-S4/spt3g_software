# Locate Python

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.12)
	find_package(Python COMPONENTS Interpreter Development REQUIRED)
else()
	find_package(PythonInterp REQUIRED)
	find_package(PythonLibs ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR} REQUIRED)
	string(REGEX REPLACE ".*libpython([0-9])\\.[0-9]+.*\\..*" "\\1" Python_VERSION_MAJOR ${PYTHON_LIBRARIES})
	string(REGEX REPLACE ".*libpython[0-9]\\.([0-9]+).*\\..*" "\\1" Python_VERSION_MINOR ${PYTHON_LIBRARIES})
	set(Python_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH})
	set(Python_LIBRARIES ${PYTHON_LIBRARIES})
	set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
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
message(STATUS "Found Boost: ${Boost_INCLUDE_DIRS} (found version \"${Boost_VERSION}\")")
