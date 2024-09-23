# Locate Python

find_package(Python COMPONENTS Interpreter Development REQUIRED)

# suppress configuration warnings in newer cmake / boost versions
if(NOT DEFINED CMAKE_FIND_PACKAGE_PREFER_CONFIG)
	set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)
endif()

# Boost bits we need
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

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
