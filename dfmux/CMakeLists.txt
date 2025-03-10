find_package(NetCDF MODULE)
find_package(HDF5)
if (NetCDF_FOUND AND HDF5_FOUND)
	set(DFMUX_LIB_EXTRA_SRC src/NetCDFDump.cxx)
else()
	set(DFMUX_LIB_EXTRA_SRC "")
endif()

add_spt3g_library(dfmux SHARED
	src/DfMuxBuilder.cxx src/DfMuxCollector.cxx src/DfMuxSample.cxx
	src/LegacyDfMuxCollector.cxx src/HardwareMap.cxx src/DfMuxCollator.cxx
	src/Housekeeping.cxx
	${DFMUX_LIB_EXTRA_SRC}
)
target_link_libraries(dfmux PUBLIC core)

add_spt3g_module(dfmux src/python.cxx)

if (NetCDF_FOUND AND HDF5_FOUND)
	target_link_libraries(dfmux PRIVATE NetCDF::NetCDF)
	add_spt3g_program(bin/ledgerman.py ledgerman)
endif()

link_python_dir()

add_spt3g_test(imports)
