if(NOT BUILD_WHEELS)
	add_executable(cppexample cppexample.cxx)
	target_link_libraries(cppexample core)
endif()
