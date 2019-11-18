add_library(maps SHARED
  src/chealpix.c
  src/get_interp_val.c
  src/G3SkyMap.cxx
  src/HealpixSkyMap.cxx
  src/FlatSkyProjection.cxx
  src/FlatSkyMap.cxx
  src/mapdata.cxx
  src/maputils.cxx
  src/pointing.cxx
  src/python.cxx
)

target_link_libraries(maps core)

link_python_dir()

add_spt3g_test(quatangtest)
add_spt3g_test(transtest SLOWTEST)
add_spt3g_test(flatsky_maps)
add_spt3g_test(flatskymap_operators)
add_spt3g_test(healpixskymap_operators)
add_spt3g_test(weights_operators)
add_spt3g_test(healpix_maps)
add_spt3g_test(fitsio)
