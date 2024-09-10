if(NOT BUILD_WHEEL)
	add_executable(cppexample cppexample.cxx)
	target_link_libraries(cppexample core)
endif()
