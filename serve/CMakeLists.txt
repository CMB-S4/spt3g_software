
function (make_embedded_resource var_name file_name out_file)
  if (${file_name} IS_NEWER_THAN ${out_file})
    message(Generating ${out_file} " from " ${file_name})
    file (WRITE ${out_file} "// automatically generated from ${file_name}\n")
    file (READ  ${file_name} filedata HEX)
    string (REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," cfiledata ${filedata})
    file (APPEND ${out_file} "unsigned char ${var_name}[] = { ${cfiledata} };\n")
  endif()
endfunction()

make_embedded_resource(favicon favicon.ico ${CMAKE_BINARY_DIR}/serve/favicon.h)

add_executable(spt3g-json-serve-files src/spt3g-json-serve-files.cxx)

target_link_libraries(spt3g-json-serve-files core z)
target_include_directories(spt3g-json-serve-files PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
)

