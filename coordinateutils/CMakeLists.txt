if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DOPENMP_FOUND)
endif()

add_library(coordinateutils SHARED
  src/chealpix.c
  src/get_interp_val.c
  src/G3SkyMap.cxx
  src/HealpixSkyMap.cxx
  src/flatskyprojection.cxx
  src/FlatSkyMap.cxx
  src/mapdata.cxx
  src/maputils.cxx
  src/pointing.cxx
  src/python.cxx
)

target_link_libraries(coordinateutils core)
execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR}/python ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/coordinateutils)

add_spt3g_test(quatangtest)
add_spt3g_test(transtest SLOWTEST)
add_spt3g_test(flatsky_maps)
add_spt3g_test(flatskymap_operators)
add_spt3g_test(healpixskymap_operators)
add_spt3g_test(weights_operators)
add_spt3g_test(healpix_maps)
add_spt3g_test(fitsio)
