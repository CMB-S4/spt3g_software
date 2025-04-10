# Convenience macros
macro(link_python_dir)
	execute_process(COMMAND mkdir -p ${SPT3G_MODULE_DIR})
	if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/python)
		execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR}/python ${SPT3G_MODULE_DIR}/${PROJECT})
		list(APPEND SPT3G_PYTHON_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/python)
	else(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/python)
		execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR} ${SPT3G_MODULE_DIR}/${PROJECT})
		list(APPEND SPT3G_PYTHON_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
	endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/python)
	set(SPT3G_PYTHON_DIRS ${SPT3G_PYTHON_DIRS} PARENT_SCOPE)
endmacro(link_python_dir)

macro(add_spt3g_program prog_name)
	set(prog_in ${CMAKE_CURRENT_SOURCE_DIR}/${prog_name})
	set(extra_macro_args ${ARGN})
	list(LENGTH extra_macro_args num_extra_args)
	if (${num_extra_args} GREATER 0)
		list(GET extra_macro_args 0 prog_out)
	else()
		get_filename_component(prog_out ${prog_in} NAME)
	endif()
	list(APPEND SPT3G_PROGRAMS ${prog_out})
	set(SPT3G_PROGRAMS ${SPT3G_PROGRAMS} PARENT_SCOPE)
	set(prog_out ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${prog_out})
	execute_process(COMMAND ln -fsn ${prog_in} ${prog_out})
endmacro(add_spt3g_program prog_name)

macro(add_spt3g_executable prog_name)
	add_executable(${prog_name} ${ARGN})
	target_link_libraries(${prog_name} pybind11::embed)
endmacro(add_spt3g_executable prog_name)

macro(add_spt3g_library lib_name)
	add_library(${lib_name} ${ARGN})
	set_target_properties(${lib_name} PROPERTIES PREFIX "libspt3g-")
	if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/include)
		target_include_directories(${lib_name} PUBLIC
			$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
	endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/include)
	# remove old libraries to avoid linking against the wrong one
	file(REMOVE ${SPT3G_MODULE_DIR}/libspt3g-${lib_name}.so)
	file(REMOVE ${SPT3G_MODULE_DIR}/libspt3g-${lib_name}.dylib)
	list(APPEND SPT3G_LIBRARIES ${lib_name})
	set(SPT3G_LIBRARIES ${SPT3G_LIBRARIES} PARENT_SCOPE)
	set_target_properties(${lib_name} PROPERTIES BUILD_RPATH "\$ORIGIN")
	install(TARGETS ${lib_name} EXPORT ${PROJECT_NAME}Config LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
endmacro(add_spt3g_library lib_name)

macro(add_spt3g_module lib_name)
	set(mod_name "_lib${lib_name}")
	pybind11_add_module(${mod_name} MODULE WITH_SOABI ${ARGN})
	if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		# Assume Linux-style ld linker
		target_link_options(${mod_name} PUBLIC "LINKER:--no-as-needed")
	endif()
	target_link_libraries(${mod_name} PUBLIC ${lib_name})
	set_target_properties(${mod_name} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${SPT3G_MODULE_DIR})
	set_target_properties(${mod_name} PROPERTIES BUILD_RPATH "\$ORIGIN/../lib")
endmacro(add_spt3g_module lib_name)

macro(add_spt3g_test test_name)
	if(DEFINED ENV{CIBUILDWHEEL})
		# The cibuildwheel builder runs tests in a separate virtual
		# environment with its own python executable that is not known
		# ahead of time.  Detect this, and assume that the correct
		# executable is already on the PATH.
		set(TEST_PYEXEC python)
	else()
		set(TEST_PYEXEC ${Python_EXECUTABLE})
	endif()
	add_test(${PROJECT}/${test_name} ${TEST_PYEXEC} ${CMAKE_CURRENT_SOURCE_DIR}/tests/${test_name}.py)

	set(extra_macro_args ${ARGN})
	list(LENGTH extra_macro_args num_extra_args)
	if(NOT DEFINED ENV{CIBUILDWHEEL})
		# install wheel before testing
		set_tests_properties(${PROJECT}/${test_name} PROPERTIES ENVIRONMENT
			"PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH};PYTHONPATH=${CMAKE_BINARY_DIR}:$ENV{PYTHONPATH};LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib:$ENV{LD_LIBRARY_PATH}")
	endif()
	if (${num_extra_args} GREATER 0)
		list(GET extra_macro_args 0 test_labels)
		set_tests_properties(${PROJECT}/${test_name} PROPERTIES LABELS ${test_labels})
	endif ()
endmacro(add_spt3g_test test_name)

if(DEFINED ENV{CIBUILDWHEEL})
	# Can't compile binaries against libpython in CI
	macro(add_spt3g_test_program test_name)
	endmacro(add_spt3g_test_program test_name)
else()
macro(add_spt3g_test_program test_name)
	cmake_parse_arguments("ADD_TEST_PROGRAM"
	                      "" # options
	                      "" # one value keywords
	                      "SOURCE_FILES;TEST_LABELS;USE_PROJECTS" # multi-value arguments
	                      ${ARGN}
	                      )

	add_test(NAME ${PROJECT}/${test_name}
	         WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
	         COMMAND ${PROJECT}-${test_name}
	         )

	if(ADD_TEST_PROGRAM_TEST_LABELS)
		set_tests_properties(${PROJECT}/${test_name} PROPERTIES LABELS ${ADD_TEST_PROGRAM_TEST_LABELS})
	endif(ADD_TEST_PROGRAM_TEST_LABELS)

	if(NOT EXISTS ${PROJECT_BINARY_DIR}/Spt3gTestMain.cxx)
		file(WRITE ${PROJECT_BINARY_DIR}/Spt3gTestMain.cxx "#include <G3Test.h>\nG3TEST_MAIN_IMPL\n")
	endif(NOT EXISTS ${PROJECT_BINARY_DIR}/Spt3gTestMain.cxx)

	add_spt3g_executable(${PROJECT}-${test_name}
	                     ${PROJECT_BINARY_DIR}/Spt3gTestMain.cxx
	                     ${ADD_TEST_PROGRAM_SOURCE_FILES}
	                     )
	target_include_directories(${PROJECT}-${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/cmake)

	foreach(USED_PROJECT ${ADD_TEST_PROGRAM_USE_PROJECTS})
		target_link_libraries(${PROJECT}-${test_name} ${USED_PROJECT})
	endforeach(USED_PROJECT ${ADD_TEST_PROGRAM_USE_PROJECTS})
endmacro(add_spt3g_test_program test_name)
endif()
