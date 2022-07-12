add_spt3g_library(maps SHARED
  src/chealpix.c
  src/G3SkyMap.cxx
  src/G3SkyMapMask.cxx
  src/HealpixSkyMapInfo.cxx
  src/HealpixSkyMap.cxx
  src/FlatSkyProjection.cxx
  src/FlatSkyMap.cxx
  src/MapBinner.cxx
  src/MapMockObserver.cxx
  src/MapTODMasker.cxx
  src/SingleDetectorBoresightBinner.cxx
  src/SingleDetectorMapBinner.cxx
  src/mapdata.cxx
  src/maputils.cxx
  src/pointing.cxx
  src/python.cxx
)

target_link_libraries(maps core calibration)
link_openmp(maps)

link_python_dir()

add_spt3g_test(quatangtest)
add_spt3g_test(transtest SLOWTEST)
add_spt3g_test(flatsky_maps)
add_spt3g_test(flatskymap_operators)
add_spt3g_test(healpixskymap_operators)
add_spt3g_test(weights_operators)
add_spt3g_test(healpix_maps)
add_spt3g_test(fitsio)
add_spt3g_test(map_modules)
