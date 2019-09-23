#include <pybindings.h>

#include <vector>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <G3Logging.h>
#include <G3Units.h>

#include <coordinateutils/HealpixSkyMap.h>
#include <coordinateutils/FlatSkyMap.h>

#include <coordinateutils/maputils.h>

#include <iostream>

using namespace G3Units;

void get_ra_dec_map_cpp(G3SkyMapConstPtr m, G3SkyMapPtr ra, G3SkyMapPtr dec)
{

	if (!m->IsCompatible(*ra)) {
		log_fatal("Output ra map must match coordinates and dimensions of input sky map");
	}
	if (!m->IsCompatible(*dec)) {
		log_fatal("Output ra map must match coordinates and dimensions of input sky map");
	}

	{
		// These are going to be dense maps, so just start that way
		FlatSkyMapPtr xra = boost::dynamic_pointer_cast<FlatSkyMap>(ra);
		FlatSkyMapPtr xdec = boost::dynamic_pointer_cast<FlatSkyMap>(dec);

		if (xra)
			xra->ConvertToDense();
		if (xdec)
			xdec->ConvertToDense();
	}

	{
		// These are going to be dense maps, so just start that way
		HealpixSkyMapPtr xra = boost::dynamic_pointer_cast<HealpixSkyMap>(ra);
		HealpixSkyMapPtr xdec = boost::dynamic_pointer_cast<HealpixSkyMap>(dec);

		if (xra)
			xra->ConvertToDense();
		if (xdec)
			xdec->ConvertToDense();
	}

	for (size_t i = 0; i < m->size(); i++) {
		std::vector<double> radec = m->pixel_to_angle(i);
		(*ra)[i] = radec[0];
		(*dec)[i] = radec[1];
	}

	ra->is_weighted = false;
	ra->units = G3Timestream::None;
	ra->pol_type = G3SkyMap::None;
	dec->is_weighted = false;
	dec->units = G3Timestream::None;
	dec->pol_type = G3SkyMap::None;

	// XXX return tuple directly
}


void reproj_map(G3SkyMapConstPtr in_map, G3SkyMapPtr out_map, int rebin, bool interp)
{

	if (in_map->coord_ref != out_map->coord_ref) {
		log_fatal("Input and output maps must use the same coordinates");
	}

	for (size_t i = 0; i < out_map->size(); i++) {
		double val = 0;
		if (rebin > 1) {
			std::vector<double> ra, dec;
			out_map->get_rebin_angles(i, rebin, ra, dec);
			for (size_t j = 0; j < ra.size(); j++) {
				if (interp)
					val += in_map->get_interp_value(ra[j], dec[j]);
				else
					val += (*in_map)[in_map->angle_to_pixel(ra[j], dec[j])];
			}
			val /= ra.size();
		} else {
			std::vector<double> radec = out_map->pixel_to_angle(i);
			if (interp)
				val = in_map->get_interp_value(radec[0], radec[1]);
			else
				val = (*in_map)[in_map->angle_to_pixel(radec[0], radec[1])];
		}
		if (val != 0)
			(*out_map)[i] = val;
	}

	out_map->coord_ref = in_map->coord_ref;
	out_map->is_weighted = in_map->is_weighted;
	out_map->units = in_map->units;
	out_map->pol_type = in_map->pol_type;
}


namespace bp = boost::python;
void maputils_pybindings(void){
	bp::def("get_ra_dec_map_cpp", get_ra_dec_map_cpp);
	bp::def("reproj_map", reproj_map,
		(bp::arg("in_map"), bp::arg("out_map"), bp::arg("rebin")=1, bp::arg("interp")=false),
		"Takes the data in in_map and reprojects it onto out_map.  out_map can\n"
		"have a different projection, size, resolution, etc.  Optionally account\n"
		"for sub-pixel structure by setting rebin > 1 and/or enable bilinear\n"
		"interpolation of values from the input map by setting interp=True");
}
