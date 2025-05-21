set(CMAKE_SOURCE_DIR ${CMAKE_ARGV3})
set(CMAKE_BINARY_DIR ${CMAKE_ARGV4})
set(Python_EXECUTABLE ${CMAKE_ARGV5})
file(GLOB cmake_projects RELATIVE ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/*/CMakeLists.txt)
foreach(d ${cmake_projects})
	get_filename_component(proj ${d} PATH)
	get_filename_component(pname ${proj} NAME_WE)
	# Skip any excluded projects
	if(EXISTS "${CMAKE_SOURCE_DIR}/${pname}/.nodocs")
		file(REMOVE "${CMAKE_SOURCE_DIR}/doc/moddoc_${pname}.rst")
		continue()
	endif()
	# Copy header if exists, or create empty rst
	if(EXISTS "${CMAKE_SOURCE_DIR}/${pname}/README.rst")
		execute_process(COMMAND ln -fsn ${CMAKE_SOURCE_DIR}/${pname}/README.rst ${CMAKE_SOURCE_DIR}/doc/intro_${pname}.rst)
		file(WRITE "${CMAKE_SOURCE_DIR}/doc/moddoc_${pname}.rst" ".. include:: intro_${pname}.rst\n\n")
	else()
		# Add an RST title. Cmake doesn't have the ability to generate
		# a string of N dashes, so generate a random string in which
		# the only character it is allowed to pick is a -
		string(LENGTH ${pname} pname_len)
		string(RANDOM LENGTH ${pname_len} ALPHABET - pname_header)
		file(WRITE "${CMAKE_SOURCE_DIR}/doc/moddoc_${pname}.rst" "${pname_header}\n${pname}\n${pname_header}\n\n")
	endif()
	execute_process(COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}:$ENV{PYTHONPATH} LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib:$ENV{LD_LIBRARY_PATH}
		${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/core/bin/spt3g-inspect "spt3g.${pname}" OUTPUT_VARIABLE MOD_DOC)
	file(APPEND "${CMAKE_SOURCE_DIR}/doc/moddoc_${pname}.rst" "${MOD_DOC}")
endforeach()
