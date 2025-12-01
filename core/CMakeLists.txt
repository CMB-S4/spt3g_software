# OS X has a broken implementatin of pthreads.
if(APPLE)
	set(CORE_EXTRA_SRCS src/ApplePthreadBarrier.cxx)
else(APPLE)
	set(CORE_EXTRA_SRCS "")
endif(APPLE)

if(NOT DEFINED WITHOUT_SUPERTIMESTREAM)
	set(WITHOUT_SUPERTIMESTREAM FALSE CACHE BOOL "Build G3SuperTimestream class")
endif()
if(NOT WITHOUT_SUPERTIMESTREAM)
	set(CORE_EXTRA_SRCS ${CORE_EXTRA_SRCS} src/G3SuperTimestream.cxx)
endif()

add_spt3g_library(core SHARED
	src/G3Logging.cxx src/G3PrintfLogger.cxx src/G3SyslogLogger.cxx
	src/crc32.c src/G3Frame.cxx src/G3Data.cxx src/G3TimeStamp.cxx
	src/G3Pipeline.cxx src/G3PipelineInfo.cxx
	src/int_storage.cxx src/G3Vector.cxx src/G3Map.cxx
	src/dataio.cxx src/compression.cxx
	src/G3Reader.cxx src/G3Writer.cxx src/G3MultiFileWriter.cxx
	src/G3NetworkSender.cxx src/G3InfiniteSource.cxx
	src/G3EventBuilder.cxx src/G3TriggeredBuilder.cxx
	src/G3Timestream.cxx src/G3Timesample.cxx src/G3Quat.cxx
	src/G3Units.cxx ${CORE_EXTRA_SRCS} src/pybindings.cxx
)


# json output option
set(ENABLE_JSON_OUTPUT "FALSE" CACHE BOOL "Enable JSON output of frames (adds to compile time and binary size)")

if (${ENABLE_JSON_OUTPUT})
	target_compile_definitions(core PUBLIC -DSPT3G_ENABLE_JSON_OUTPUT)
endif()


# Link dependencies
target_link_libraries(core PUBLIC spt3g)

find_package(Filesystem REQUIRED)
target_link_libraries(core PRIVATE std::filesystem)

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
	target_link_libraries(core PRIVATE FLAC::FLAC)
endif()

if(NOT WITHOUT_SUPERTIMESTREAM)
	find_package(OpenMP QUIET)
	if(OpenMP_FOUND)
		target_link_libraries(core PRIVATE OpenMP::OpenMP_CXX)
	endif()
endif()

# Link against Z library
if(NOT DEFINED WITH_ZLIB)
	set(WITH_ZLIB TRUE CACHE BOOL "Enable gzip fie compression")
endif()
if(WITH_ZLIB)
	find_package(ZLIB)
endif()
if(ZLIB_FOUND)
	target_compile_definitions(core PRIVATE -DZLIB_FOUND)
	target_link_libraries(core PRIVATE ZLIB::ZLIB)
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
	target_link_libraries(core PRIVATE BZip2::BZip2)
endif()

if(NOT DEFINED WITH_LZMA)
	set(WITH_LZMA TRUE CACHE BOOL "Enable LZMA file compression")
endif()
if(WITH_LZMA)
	find_package(LibLZMA)
endif()
if(LIBLZMA_FOUND)
	target_compile_definitions(core PRIVATE -DLZMA_FOUND)
	target_link_libraries(core PRIVATE LibLZMA::LibLZMA)
endif()

link_python_dir()

add_spt3g_program(bin/spt3g-dump)
add_spt3g_program(bin/spt3g-jsonify)
add_spt3g_program(bin/spt3g-verify)
add_spt3g_program(bin/spt3g-inspect)
add_spt3g_program(bin/gen-analysis-doc)

#Tests
add_spt3g_test(imports)
add_spt3g_test(copycons)
add_spt3g_test(maps)
add_spt3g_test(framepickle)
add_spt3g_test(pipeline)
add_spt3g_test(pipeline_module)
add_spt3g_test(multiproc)
add_spt3g_test(frames)
add_spt3g_test(frame_types)
add_spt3g_test(frameobjects)
add_spt3g_test(datatypes)
if(FLAC_FOUND)
	add_spt3g_test(ts_nanencoding)
	add_spt3g_test(ts_flac_compression)
endif()
add_spt3g_test(cuts)
add_spt3g_test(fileio)
add_spt3g_test(multifileio)
add_spt3g_test(splitfileio)
if(ZLIB_FOUND)
	add_spt3g_test(compressedfileio)
endif()
if(BZIP2_FOUND)
	add_spt3g_test(bz2fileio)
endif()
if(LIBLZMA_FOUND)
	add_spt3g_test(xzfileio)
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
add_spt3g_test(non_contiguous_vectors)
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
                         ${CMAKE_CURRENT_SOURCE_DIR}/tests/OrderedMapTest.cxx
                       USE_PROJECTS core)
