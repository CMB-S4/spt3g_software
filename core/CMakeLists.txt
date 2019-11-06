# OS X has a broken implementatin of pthreads.
if(APPLE)
	set(CORE_EXTRA_SRCS src/ApplePthreadBarrier.cxx)
else(APPLE)
	set(CORE_EXTRA_SRCS "")
endif(APPLE)

add_spt3g_library(core SHARED
	src/G3EventBuilder.cxx src/G3Frame.cxx src/G3TimeStamp.cxx
	src/G3Pipeline.cxx src/G3Writer.cxx src/G3Reader.cxx src/python.cxx
	src/G3InfiniteSource.cxx src/G3Logging.cxx src/G3PrintfLogger.cxx
	src/G3Data.cxx src/G3Vector.cxx src/G3Map.cxx src/G3Timestream.cxx
        src/G3Timesample.cxx
	src/G3TriggeredBuilder.cxx src/G3MultiFileWriter.cxx src/dataio.cxx
	src/crc32.c ${CORE_EXTRA_SRCS}
	src/G3NetworkSender.cxx src/G3SyslogLogger.cxx
	src/G3PipelineInfo.cxx src/G3Quat.cxx
)

# Ignore known warnings
target_compile_options(core PUBLIC -Wno-unknown-warning-option -Wno-unused -Wno-unused-result -Wno-deprecated-register -Wno-sign-compare -Wno-constant-conversion -Wno-self-assign-overloaded -Wno-deprecated-declarations)

# Fix bugs in GCC 4.4's strict aliasing code by turning it off
target_compile_options(core PUBLIC -fno-strict-aliasing)

# Make core includes available without directories
target_include_directories(core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/core)

# Link against FLAC library
find_package(FLAC)
if(FLAC_FOUND)
	target_compile_definitions(core PRIVATE -DG3_HAS_FLAC)
	target_include_directories(core PRIVATE ${FLAC_INCLUDE_DIR})
	target_link_libraries(core ${FLAC_LIBRARIES})
endif()

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
if(Boost_VERSION EQUAL 104700 OR Boost_VERSION GREATER 104700)
	set(Boost_GEOMETRY_AVAIL ON)
	message(STATUS "Boost Geometry Available")
endif()

target_include_directories(core PUBLIC ${Boost_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH})
target_link_libraries(core pthread ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

# Work around yet more bugs in GCC 4.4, this time with C++ 11 support
target_compile_definitions(core PUBLIC -DBOOST_NO_CXX11_SMART_PTR=1)

link_python_dir()

execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR}/bin/spt3g-dump ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR}/bin/spt3g-verify ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR}/bin/spt3g-inspect ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR}/bin/gen-analysis-doc ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

#Tests
add_spt3g_test(copycons)
add_spt3g_test(framepickle)
add_spt3g_test(pipeline)
add_spt3g_test(pipeline_module)
add_spt3g_test(multiproc)
add_spt3g_test(frames)
add_spt3g_test(frame_types)
add_spt3g_test(datatypes)
if(FLAC_FOUND)
	add_spt3g_test(ts_nanencoding)
endif()
add_spt3g_test(cuts)
add_spt3g_test(fileio)
add_spt3g_test(multifileio)
add_spt3g_test(splitfileio)
add_spt3g_test(compressedfileio)
add_spt3g_test(portability)
add_spt3g_test(ts_bufferprotocol)
add_spt3g_test(timestream_slicing)
add_spt3g_test(timestream_times)
add_spt3g_test(deduplicate)
add_spt3g_test(timestrings)
add_spt3g_test(testcomplexconv)
add_spt3g_test(droporphanmetadata)
add_spt3g_test(networktest SLOWTEST)
add_spt3g_test(endprocessing)
add_spt3g_test(lazyreader)
add_spt3g_test(returnvals)
add_spt3g_test(timestream_2d)
add_spt3g_test(cached_frame)
add_spt3g_test(pipelineinfo)
add_spt3g_test(quaternions)
add_spt3g_test(timesample)
