if(NOT DEFINED ENV{CIBUILDWHEEL})
	# Can't compile binaries against libpython in CI
	add_spt3g_executable(cppexample cppexample.cxx)
	target_link_libraries(cppexample core)
endif()
