#include <pybindings.h>

#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <G3Logging.h>
#include <G3Units.h>

#include <maps/maputils.h>

#include <iostream>

#define COS cos
#define SIN sin
#define ATAN2 atan2

using namespace G3Units;

G3SkyMapPtr GetMaskMap(G3SkyMapConstPtr m)
{
	G3SkyMapPtr mask = m->Clone(false);

	for (size_t i = 0; i < m->size(); i++) {
		if (m->at(i) == 0)
			continue;
		(*mask)[i] = 1.0;
	}

	return mask;
}

void RemoveWeights(G3SkyMapPtr T, G3SkyMapPtr Q, G3SkyMapPtr U, G3SkyMapWeightsConstPtr W,
    bool zero_nans)
{
	bool pol = W->IsPolarized();

	g3_assert(T->weighted);
	g3_assert(W->IsCongruent());
	g3_assert(T->IsCompatible(*(W->TT)));

	if (pol) {
		g3_assert(!!Q && !!U);
		g3_assert(T->IsCompatible(*Q));
		g3_assert(T->IsCompatible(*U));
		g3_assert(Q->weighted);
		g3_assert(U->weighted);
	}

	if (!zero_nans) {
		T->ConvertToDense();
		if (pol) {
			Q->ConvertToDense();
			U->ConvertToDense();
			for (size_t pix = 0; pix < T->size(); pix++) {
				StokesVector v((*T)[pix], (*Q)[pix], (*U)[pix]);
				v /= W->at(pix);
			}
		} else {
			(*T) /= *(W->TT);
		}
	} else {
		for (size_t pix = 0; pix < W->size(); pix++) {
			double v = T->at(pix);
			MuellerMatrix m = W->at(pix);
			bool empty;

			// skip empty pixels
			if (!pol) {
				empty = (m.tt == 0);
				if (empty && v == 0)
					continue;
			} else {
				empty = (m.Det() < 1e-12);
				if (empty && v == 0 && Q->at(pix) == 0 && U->at(pix) == 0)
					continue;
			}

			// set bad pixels to 0
			if (empty) {
				(*T)[pix] = 0;
				if (pol) {
					(*Q)[pix] = 0;
					(*U)[pix] = 0;
				}
				continue;
			}

			// remove weights
			if (pol) {
				StokesVector v((*T)[pix], (*Q)[pix], (*U)[pix]);
				v /= m;
			} else {
				(*T)[pix] /= (*(W->TT))[pix];
			}
		}
	}

	T->weighted = false;
	if (pol) {
		Q->weighted = false;
		U->weighted = false;
	}
}

void ApplyWeights(G3SkyMapPtr T, G3SkyMapPtr Q, G3SkyMapPtr U, G3SkyMapWeightsConstPtr W)
{
	bool pol = W->IsPolarized();

	g3_assert(!T->weighted);
	g3_assert(W->IsCongruent());
	g3_assert(T->IsCompatible(*(W->TT)));

	if (pol) {
		g3_assert(!!Q && !!U);
		g3_assert(T->IsCompatible(*Q));
		g3_assert(T->IsCompatible(*U));
		g3_assert(!Q->weighted);
		g3_assert(!U->weighted);
	}

	if (pol) {
		for (size_t pix = 0; pix < T->size(); pix++) {
			if (T->at(pix) == 0 && Q->at(pix) == 0 && U->at(pix) == 0)
				continue;
			StokesVector v((*T)[pix], (*Q)[pix], (*U)[pix]);
			v = W->at(pix) * v;
		}
	} else {
		(*T) *= *(W->TT);
	}

	T->weighted = true;
	if (pol) {
		Q->weighted = true;
		U->weighted = true;
	}
}

boost::python::tuple GetRaDecMap(G3SkyMapConstPtr m)
{

	G3SkyMapPtr ra = m->Clone(false);
	G3SkyMapPtr dec = m->Clone(false);

	// These are going to be dense maps, so just start that way
	ra->ConvertToDense();
	dec->ConvertToDense();

	for (size_t i = 0; i < m->size(); i++) {
		std::vector<double> radec = m->PixelToAngle(i);
		(*ra)[i] = radec[0];
		(*dec)[i] = radec[1];
	}

	ra->weighted = false;
	ra->units = G3Timestream::None;
	ra->pol_type = G3SkyMap::None;
	dec->weighted = false;
	dec->units = G3Timestream::None;
	dec->pol_type = G3SkyMap::None;

	return boost::python::make_tuple(ra, dec);
}


void FlattenPol(FlatSkyMapPtr Q, FlatSkyMapPtr U, double h, bool invert)
{
	g3_assert(Q->IsCompatible(*U));
	g3_assert(Q->IsPolFlat() == U->IsPolFlat());

	if (Q->IsPolFlat() && !invert)
		return;
	if (!Q->IsPolFlat() && invert)
		return;

	for (auto i : *Q) {
		double q = i.second;
		double u = U->at(i.first);
		if (q == 0 && u == 0)
			continue;

		std::vector<double> grad = Q->PixelToAngleGrad(i.first, h);
		double rot = ATAN2(-grad[0], grad[1]) + ATAN2(-grad[3], -grad[2]);
		if (invert)
			rot *= -1.0;
		double cr = COS(rot);
		double sr = SIN(rot);

		(*Q)[i.first] = cr * q - sr * u;
		(*U)[i.first] = sr * q + cr * u;
	}

	Q->SetFlatPol(!invert);
	U->SetFlatPol(!invert);
}


void ReprojMap(G3SkyMapConstPtr in_map, G3SkyMapPtr out_map, int rebin, bool interp)
{

	if (in_map->coord_ref != out_map->coord_ref) {
		log_fatal("Input and output maps must use the same coordinates");
	}

	for (size_t i = 0; i < out_map->size(); i++) {
		double val = 0;
		if (rebin > 1) {
			std::vector<double> ra, dec;
			out_map->GetRebinAngles(i, rebin, ra, dec);
			for (size_t j = 0; j < ra.size(); j++) {
				if (interp)
					val += in_map->GetInterpValue(ra[j], dec[j]);
				else
					val += in_map->at(in_map->AngleToPixel(ra[j], dec[j]));
			}
			val /= ra.size();
		} else {
			std::vector<double> radec = out_map->PixelToAngle(i);
			if (interp)
				val = in_map->GetInterpValue(radec[0], radec[1]);
			else
				val = in_map->at(in_map->AngleToPixel(radec[0], radec[1]));
		}
		if (val != 0)
			(*out_map)[i] = val;
	}

	out_map->coord_ref = in_map->coord_ref;
	out_map->weighted = in_map->weighted;
	out_map->units = in_map->units;
	out_map->pol_type = in_map->pol_type;
}


namespace bp = boost::python;
void maputils_pybindings(void){
	bp::def("get_mask_map", GetMaskMap, (bp::arg("map_in")),
		"Returns a map that is 1 where the input map is nonzero.");

	bp::def("remove_weights", RemoveWeights,
		(bp::arg("T"), bp::arg("Q"), bp::arg("U"), bp::arg("W"), bp::arg("zero_nans")=false),
		"Remove weights from polarized or unpolarized maps. "
		"If zero_nans is true, empty pixels are skipped, and pixels "
		"with zero weight are set to 0 instead of nan.");

	bp::def("apply_weights", ApplyWeights,
		(bp::arg("T"), bp::arg("Q"), bp::arg("U"), bp::arg("W")),
		"Remove weights from polarized or unpolarized maps.");

	bp::def("get_ra_dec_map", GetRaDecMap, (bp::arg("map_in")),
		"Returns maps of the ra and dec angles for each pixel in the input map");

	bp::def("flatten_pol", FlattenPol,
		(bp::arg("Q"), bp::arg("U"), bp::arg("h")=0.001, bp::arg("invert")=false),
		"For maps defined on the sphere the direction of the polarization angle is "
		"is defined relative to the direction of North.  When making maps we follow "
		"this definition.\n\nFor any flat sky estimators, the polarization angle is "
		"defined relative to the vertical axis.  For some map projections the "
		"direction of north is not the same as the vertical axis.  This function "
		"applies a rotation to the Q and U values to switch the curved sky Q/U "
		"definition to the flat sky Q/U definition.\n\nIf for whatever reason you "
		"want to reverse the process set the invert argument to True.");

	bp::def("reproj_map", ReprojMap,
		(bp::arg("in_map"), bp::arg("out_map"), bp::arg("rebin")=1, bp::arg("interp")=false),
		"Takes the data in in_map and reprojects it onto out_map.  out_map can\n"
		"have a different projection, size, resolution, etc.  Optionally account\n"
		"for sub-pixel structure by setting rebin > 1 and/or enable bilinear\n"
		"interpolation of values from the input map by setting interp=True");
}
