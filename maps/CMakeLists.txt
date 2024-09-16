add_spt3g_library(maps
  src/chealpix.c
  src/G3SkyMap.cxx
  src/G3SkyMapMask.cxx
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
  src/mapdata.cxx
  src/maputils.cxx
  src/pointing.cxx
  src/python.cxx
)

if(BUILD_WHEEL)
	target_include_directories(maps PUBLIC ${CMAKE_SOURCE_DIR}/calibration/include)
else()
	target_link_libraries(maps PUBLIC core calibration)
endif()
link_openmp(maps)

link_python_dir()

add_spt3g_program(bin/spt3g-coadd-maps)

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
