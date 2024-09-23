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
set(Boost_PYTHON_VERSION ${Python_VERSION})

find_package(Boost COMPONENTS system iostreams filesystem python REQUIRED)
message(STATUS "Found Boost: ${Boost_INCLUDE_DIRS} (found version \"${Boost_VERSION}\")")
