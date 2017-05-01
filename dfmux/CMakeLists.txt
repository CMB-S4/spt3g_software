find_package(NetCDF)
if (NETCDF_FOUND)
	set(DFMUX_LIB_EXTRA_SRC src/NetCDFDump.cxx)
	set(DFMUX_LIB_EXTRA_LIB ${NETCDF_LIBRARIES})
else()
	set(DFMUX_LIB_EXTRA_SRC "")
	set(DFMUX_LIB_EXTRA_LIB "")
endif()

add_library(dfmux SHARED
	src/DfMuxBuilder.cxx src/DfMuxCollector.cxx src/DfMuxSample.cxx
	src/LegacyDfMuxCollector.cxx src/HardwareMap.cxx src/DfMuxCollator.cxx
	src/Housekeeping.cxx src/python.cxx
	${DFMUX_LIB_EXTRA_SRC}
)
target_link_libraries(dfmux core ${DFMUX_LIB_EXTRA_LIB})

if (NETCDF_FOUND)
	execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR}/bin/ledgerman.py ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ledgerman)
endif()

link_python_dir()
