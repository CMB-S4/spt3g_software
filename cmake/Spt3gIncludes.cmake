# Convenience macros
macro(link_python_dir)
	execute_process(COMMAND mkdir -p ${SPT3G_LIBRARY_DIR})
	if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/python)
		execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR}/python ${SPT3G_LIBRARY_DIR}/${PROJECT})
		list(APPEND SPT3G_PYTHON_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/python)
	else(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/python)
		execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR} ${SPT3G_LIBRARY_DIR}/${PROJECT})
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

macro(add_spt3g_library lib_name)
	add_library(${lib_name} ${ARGN})
	if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/include)
		target_include_directories(${lib_name} PUBLIC
			$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
	endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/include)
	list(APPEND SPT3G_LIBRARIES ${lib_name})
	set(SPT3G_LIBRARIES ${SPT3G_LIBRARIES} PARENT_SCOPE)
	install(TARGETS ${lib_name} EXPORT ${PROJECT_NAME}Config LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
	install(TARGETS ${lib_name} DESTINATION ${PYTHON_MODULE_DIR}/spt3g)
endmacro(add_spt3g_library lib_name)

macro(add_spt3g_test test_name)
	add_test(${PROJECT}/${test_name} ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/${test_name}.py)

	set(extra_macro_args ${ARGN})
	list(LENGTH extra_macro_args num_extra_args)
	set_tests_properties(${PROJECT}/${test_name} PROPERTIES ENVIRONMENT
		"PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH};PYTHONPATH=${CMAKE_BINARY_DIR}:$ENV{PYTHONPATH};LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/spt3g:$ENV{LD_LIBRARY_PATH}")
	if (${num_extra_args} GREATER 0)
		list(GET extra_macro_args 0 test_labels)
		set_tests_properties(${PROJECT}/${test_name} PROPERTIES LABELS ${test_labels})
	endif ()
endmacro(add_spt3g_test test_name)

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
	
	add_executable(${PROJECT}-${test_name}
	               ${PROJECT_BINARY_DIR}/Spt3gTestMain.cxx
	               ${ADD_TEST_PROGRAM_SOURCE_FILES}
	               )
	target_include_directories(${PROJECT}-${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/cmake)
	
	foreach(USED_PROJECT ${ADD_TEST_PROGRAM_USE_PROJECTS})
		target_link_libraries(${PROJECT}-${test_name} ${USED_PROJECT})
	endforeach(USED_PROJECT ${ADD_TEST_PROGRAM_USE_PROJECTS})
endmacro(add_spt3g_test_program test_name)
