if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DOPENMP_FOUND)
endif()

add_library(coordinateutils SHARED
  src/chealpix.c
  src/wcscon.c
  src/hget.c
  src/iget.c
  src/get_interp_val.c
  src/polarizationutils.cpp
  src/coordinateutils.cpp
  src/CutskyHealpixMap.cxx
  src/FlatSkyMap.cxx
  src/python.cpp
)

target_link_libraries(coordinateutils core)
execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR}/python ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/coordinateutils)

#Tests
add_spt3g_test(coord_change_fixture SLOWTEST)
add_spt3g_test(polarization_transforms_fixture SLOWTEST)
