# Find the numpy include path
execute_process(COMMAND ${Python_EXECUTABLE} -c
  "import numpy; print(numpy.get_include())"
  OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

add_spt3g_library(so SHARED
  src/G3SuperTimestream.cxx
  src/numpy_assist.cxx
)

target_include_directories(so PUBLIC ${NUMPY_INCLUDE_DIR})
target_link_libraries(so PUBLIC core)

find_package(OpenMP)
if(OpenMP_FOUND)
  target_link_libraries(so PUBLIC OpenMP::OpenMP_CXX)
  set(OpenMP_FOUND ${OpenMP_FOUND} PARENT_SCOPE)
endif()

add_spt3g_module(so src/python.cxx)

link_python_dir()

add_spt3g_test(imports)
