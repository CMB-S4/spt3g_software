find_package(NetCDF MODULE)
if (NetCDF_FOUND)
	set(DFMUX_LIB_EXTRA_SRC src/NetCDFDump.cxx)
else()
	set(DFMUX_LIB_EXTRA_SRC "")
endif()

add_spt3g_library(dfmux SHARED
	src/DfMuxBuilder.cxx src/DfMuxCollector.cxx src/DfMuxSample.cxx
	src/LegacyDfMuxCollector.cxx src/HardwareMap.cxx src/DfMuxCollator.cxx
	src/Housekeeping.cxx src/python.cxx
	${DFMUX_LIB_EXTRA_SRC}
)
target_link_libraries(dfmux PUBLIC core)
if (NetCDF_FOUND)
	target_link_libraries(dfmux PUBLIC NetCDF::NetCDF)
endif()

if (NetCDF_FOUND)
	add_spt3g_program(bin/ledgerman.py ledgerman)
endif()

link_python_dir()
