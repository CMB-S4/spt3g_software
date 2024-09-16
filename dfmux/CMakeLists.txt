find_package(NetCDF)
if (NETCDF_FOUND)
	set(DFMUX_LIB_EXTRA_SRC src/NetCDFDump.cxx)
	set(DFMUX_LIB_EXTRA_LIB ${NETCDF_LIBRARIES})
else()
	set(DFMUX_LIB_EXTRA_SRC "")
	set(DFMUX_LIB_EXTRA_LIB "")
endif()

add_spt3g_library(dfmux
	src/DfMuxBuilder.cxx src/DfMuxCollector.cxx src/DfMuxSample.cxx
	src/LegacyDfMuxCollector.cxx src/HardwareMap.cxx src/DfMuxCollator.cxx
	src/Housekeeping.cxx src/python.cxx
	${DFMUX_LIB_EXTRA_SRC}
)
if(NOT BUILD_WHEEL)
	target_link_libraries(dfmux PUBLIC core)
endif()
if (NETCDF_FOUND)
	target_link_libraries(dfmux PUBLIC ${NETCDF_LIBRARIES})
	target_include_directories(dfmux PRIVATE ${NETCDF_INCLUDES})
endif()

if (NETCDF_FOUND)
	add_spt3g_program(bin/ledgerman.py ledgerman)
endif()

link_python_dir()
