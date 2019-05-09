#include <pybindings.h>

#include <vector>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <G3Logging.h>
#include <G3Units.h>
#include <G3SkyMap.h>

#include <coordinateutils/CutSkyHealpixMap.h>

#include <coordinateutils/maputils.h>

#include <iostream>

#ifdef OPENMP_FOUND
#include <omp.h>
#endif

using namespace G3Units;

void get_ra_dec_map_cpp(G3SkyMapConstPtr m, G3SkyMapPtr ra, G3SkyMapPtr dec)
{
	ra->EnsureAllocated();
	dec->EnsureAllocated();

	if (!m->IsCompatible(*ra)) {
		log_fatal("Output ra map must match coordinates and dimensions of input sky map");
	}
	if (!m->IsCompatible(*dec)) {
		log_fatal("Output ra map must match coordinates and dimensions of input sky map");
	}

#ifdef OPENMP_FOUND
#pragma omp parallel for
#endif
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
}


void reproj_map(G3SkyMapConstPtr in_map, G3SkyMapPtr out_map, int rebin)
{
	out_map->EnsureAllocated();
	if (in_map->coord_ref != out_map->coord_ref) {
		log_fatal("Input and output maps must use the same coordinates");
	}

#ifdef OPENMP_FOUND
#pragma omp parallel for
#endif
	for (size_t i = 0; i < out_map->size(); i++) {
		double val = 0;
		if (rebin > 1) {
			std::vector<double> ra, dec;
			out_map->get_rebin_angles(i, rebin, ra, dec);
			for (size_t j = 0; j < ra.size(); j++) {
				val += in_map->get_interp_value(ra[j], dec[j]);
			}
			val /= ra.size();
		} else {
			std::vector<double> radec = out_map->pixel_to_angle(i);
			val = in_map->get_interp_value(radec[0], radec[1]);
		}
		(*out_map)[i] = val;
	}

	out_map->coord_ref = in_map->coord_ref;
	out_map->is_weighted = in_map->is_weighted;
	out_map->units = in_map->units;
	out_map->pol_type = in_map->pol_type;
}


void reproj_fullsky_healpix_map(std::vector<double> in_map, G3SkyMapPtr out_map,
    bool nest, int rebin)
{
	//grab our out ra dec values
	out_map->EnsureAllocated();

	long nside = npix2nside(in_map.size());
	if (nside < 0) {
		log_fatal("Input map has invalid size for a healpix map");
	}
	HealpixHitPix hitpix(nside, nest, out_map->coord_ref);

#ifdef OPENMP_FOUND
#pragma omp parallel for
#endif
	for (size_t i = 0; i < out_map->size(); i++) {
		double val = 0;
		if (rebin > 1) {
			std::vector<double> ra, dec;
			out_map->get_rebin_angles(i, rebin, ra, dec);
			for (size_t j = 0; j < ra.size(); j++) {
				std::vector<long> pixels;
				std::vector<double> weights;
				hitpix.get_interp_pixels_weights(ra[j], dec[j], pixels, weights, false);
				for (size_t k = 0; k < pixels.size(); k++) {
					val += in_map[pixels[k]] * weights[k];
				}
			}
			val /= ra.size();
		} else {
			std::vector<double> radec = out_map->pixel_to_angle(i);
			std::vector<long> pixels;
			std::vector<double> weights;
			hitpix.get_interp_pixels_weights(radec[0], radec[1], pixels, weights, false);
			for (size_t j = 0; j < pixels.size(); j++) {
				val += in_map[pixels[j]] * weights[j];
			}
		}
		(*out_map)[i] = val;
	}
}


namespace bp = boost::python;
void maputils_pybindings(void){
	bp::def("get_ra_dec_map_cpp", get_ra_dec_map_cpp);
	bp::def("reproj_map", reproj_map,
		(bp::arg("in_map"), bp::arg("out_map"), bp::arg("rebin")=1),
		"Takes the data in in_map and reprojects it onto out_map.  out_map can\n"
		"have a different projection, size, resolution, etc.  Optionally account\n"
		"for sub-pixel structure in the interpolation by setting rebin > 1.");
	bp::def("reproj_fullsky_healpix_map", reproj_fullsky_healpix_map,
		( bp::arg("in_map"), bp::arg("out_map"), bp::arg("nest")=false, bp::arg("rebin")=1),
		"Takes the data in in_map (a full sky healpix map stored as a simple array)\n"
		"and reprojects it onto out_map.  out_map can be any G3SkyMap instance.\n"
		"Optionally account for sub-pixel structure in the interpolation by setting\n"
		"rebin > 1");
}
