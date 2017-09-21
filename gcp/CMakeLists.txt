add_library(gcp SHARED
	src/ARCFileReader.cxx src/python.cxx
	src/ACUStatus.cxx src/TrackerStatus.cxx src/TrackerPointing.cxx
	src/GCPMuxDataDecoder.cxx src/GCPLogger.cxx
)
target_link_libraries(gcp core dfmux)
link_python_dir()
