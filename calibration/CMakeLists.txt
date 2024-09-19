add_spt3g_library(calibration SHARED
	src/BoloProperties.cxx src/PointingProperties.cxx
)
target_link_libraries(calibration core)
link_python_dir()

add_spt3g_module(calibration src/python.cxx)

add_spt3g_test(can_modify_bolo_props_in_map SLOWTEST)
add_spt3g_test(band_format)
