set(CMAKE_SOURCE_DIR ${CMAKE_ARGV3})
set(CMAKE_BINARY_DIR ${CMAKE_ARGV4})
file(GLOB cmake_projects RELATIVE ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/*/CMakeLists.txt)
foreach(d ${cmake_projects})
	get_filename_component(proj ${d} PATH)
	get_filename_component(pname ${proj} NAME_WE)
	# Copy header if exists, or create empty rst
	if(EXISTS "${CMAKE_SOURCE_DIR}/${pname}/README.rst")
		file(READ "${CMAKE_SOURCE_DIR}/${pname}/README.rst" MOD_HEADER)
		file(WRITE "${CMAKE_SOURCE_DIR}/doc/moddoc_${pname}.rst" "${MOD_HEADER}")
	else()
		# Add an RST title. Cmake doesn't have the ability to generate
		# a string of N dashes, so generate a random string in which
		# the only character it is allowed to pick is a -
		string(LENGTH ${pname} pname_len)
		string(RANDOM LENGTH ${pname_len} ALPHABET - pname_header)
		file(WRITE "${CMAKE_SOURCE_DIR}/doc/moddoc_${pname}.rst" "${pname_header}\n${pname}\n${pname_header}\n\n")
	endif()
	execute_process(COMMAND spt3g-inspect "spt3g.${pname}" OUTPUT_VARIABLE MOD_DOC)
	file(APPEND "${CMAKE_SOURCE_DIR}/doc/moddoc_${pname}.rst" "${MOD_DOC}")
endforeach()
