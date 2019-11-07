add_spt3g_library(coordinateutils SHARED
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

set(OpenMP_FIND_QUIETLY TRUE)
find_package(OpenMP)
if(OPENMP_FOUND)
	target_compile_definitions(coordinateutils PUBLIC -DOPENMP_FOUND)
	target_compile_options(coordinateutils PUBLIC ${OpenMP_CXX_FLAGS})
endif()

target_link_libraries(coordinateutils core)

link_python_dir()

add_spt3g_test(quatangtest)
add_spt3g_test(transtest SLOWTEST)
add_spt3g_test(flatsky_maps)
add_spt3g_test(flatskymap_operators)
add_spt3g_test(healpixskymap_operators)
add_spt3g_test(weights_operators)
add_spt3g_test(healpix_maps)
add_spt3g_test(fitsio)
