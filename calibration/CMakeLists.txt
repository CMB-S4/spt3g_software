add_spt3g_library(calibration SHARED
	src/BoloProperties.cxx src/PointingProperties.cxx src/python.cxx
)
target_link_libraries(calibration PUBLIC core)
link_python_dir()

add_spt3g_test(can_modify_bolo_props_in_map SLOWTEST)
add_spt3g_test(band_format)
