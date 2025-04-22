add_spt3g_library(maps SHARED
  src/chealpix.c
  src/mapdata.cxx
  src/G3SkyMapMask.cxx
  src/G3SkyMap.cxx
  src/HealpixSkyMapInfo.cxx
  src/HealpixSkyMap.cxx
  src/FlatSkyProjection.cxx
  src/FlatSkyMap.cxx
  src/HitsBinner.cxx
  src/MapBinner.cxx
  src/MapMockObserver.cxx
  src/MapTODMasker.cxx
  src/MapTODPointing.cxx
  src/SingleDetectorBoresightBinner.cxx
  src/SingleDetectorMapBinner.cxx
  src/maputils.cxx
  src/pointing.cxx
)

target_link_libraries(maps PUBLIC core calibration)

find_package(OpenMP)
if(OpenMP_FOUND)
  target_link_libraries(maps PRIVATE OpenMP::OpenMP_CXX)
endif()

add_spt3g_module(maps src/python.cxx)

link_python_dir()

add_spt3g_program(bin/spt3g-coadd-maps)

add_spt3g_test(imports)
add_spt3g_test(quatangtest)
add_spt3g_test(transtest SLOWTEST)
add_spt3g_test(flatsky_maps)
add_spt3g_test(flatskymap_operators)
add_spt3g_test(healpixskymap_operators)
add_spt3g_test(weights_operators)
add_spt3g_test(mask_operators)
add_spt3g_test(healpix_maps)
add_spt3g_test(fitsio)
add_spt3g_test(map_modules)
add_spt3g_test(mask_serialization_test)
