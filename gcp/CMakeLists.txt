add_spt3g_library(gcp
	src/ARCFileReader.cxx src/python.cxx
	src/ACUStatus.cxx src/TrackerStatus.cxx src/TrackerPointing.cxx
	src/GCPMuxDataDecoder.cxx src/GCPLogger.cxx
)
target_link_libraries(gcp PUBLIC core dfmux)

link_python_dir()

add_spt3g_program(bin/spt3g-arc-dump)
add_spt3g_program(bin/spt3g-arc-verify)
