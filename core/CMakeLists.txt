# OS X has a broken implementatin of pthreads.
if(APPLE)
	set(CORE_EXTRA_SRCS src/ApplePthreadBarrier.cxx)
else(APPLE)
	set(CORE_EXTRA_SRCS "")
endif(APPLE)

add_spt3g_library(core SHARED
	src/G3EventBuilder.cxx src/G3Frame.cxx src/G3TimeStamp.cxx
	src/G3Pipeline.cxx src/G3Writer.cxx src/G3Reader.cxx
	src/G3InfiniteSource.cxx src/G3Logging.cxx src/G3PrintfLogger.cxx
	src/G3Data.cxx src/G3Vector.cxx src/G3Map.cxx src/G3Timestream.cxx
        src/G3Timesample.cxx
	src/G3TriggeredBuilder.cxx src/G3MultiFileWriter.cxx src/dataio.cxx
	src/crc32.c ${CORE_EXTRA_SRCS}
	src/G3NetworkSender.cxx src/G3SyslogLogger.cxx
	src/G3PipelineInfo.cxx src/G3Quat.cxx src/G3Units.cxx
	src/int_storage.cxx src/pybindings.cxx
)

# Link dependencies
target_link_libraries(core PUBLIC spt3g)

find_package(Filesystem REQUIRED)
target_link_libraries(core PUBLIC std::filesystem)

# add python bindings
add_spt3g_module(core src/python.cxx)

# Make core includes available without directories
target_include_directories(core PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/core>
	$<INSTALL_INTERFACE:${SPT3G_INCLUDE_INSTALL_DIR}/core>
)

# Link against FLAC library
find_package(FLAC MODULE)
if(FLAC_FOUND)
	target_compile_definitions(core PRIVATE -DG3_HAS_FLAC)
	target_link_libraries(core PUBLIC FLAC::FLAC)
	set(FLAC_FOUND ${FLAC_FOUND} PARENT_SCOPE)
endif()

# Link against BZIP2 library
if(NOT DEFINED WITH_BZIP2)
	set(WITH_BZIP2 TRUE CACHE BOOL "Enable bzip2 file compression")
endif()
if(WITH_BZIP2)
	find_package(BZip2)
endif()
if(BZIP2_FOUND)
	target_compile_definitions(core PRIVATE -DBZIP2_FOUND)
endif()

link_python_dir()

add_spt3g_program(bin/spt3g-dump)
add_spt3g_program(bin/spt3g-verify)
add_spt3g_program(bin/spt3g-inspect)
add_spt3g_program(bin/gen-analysis-doc)

#Tests
add_spt3g_test(imports)
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
if(BZIP2_FOUND)
	add_spt3g_test(bz2fileio)
endif()
add_spt3g_test(portability)
add_spt3g_test(vecint)
add_spt3g_test(ts_bufferprotocol)
add_spt3g_test(timestream_slicing)
add_spt3g_test(timestream_times)
add_spt3g_test(timestream_ops)
add_spt3g_test(deduplicate)
add_spt3g_test(timestrings)
add_spt3g_test(testcomplexconv)
add_spt3g_test(droporphanmetadata)
add_spt3g_test(networktest SLOWTEST)
add_spt3g_test(non_double_timestreams)
add_spt3g_test(endprocessing)
add_spt3g_test(lazyreader)
add_spt3g_test(returnvals)
add_spt3g_test(timestream_2d)
add_spt3g_test(cached_frame)
add_spt3g_test(pipelineinfo)
add_spt3g_test(quaternions)
add_spt3g_test(timesample)

add_spt3g_test_program(test
                       SOURCE_FILES
                         ${CMAKE_CURRENT_SOURCE_DIR}/tests/G3TimestreamTest.cxx
                         ${CMAKE_CURRENT_SOURCE_DIR}/tests/G3TimestreamMapTest.cxx
                       USE_PROJECTS core)
