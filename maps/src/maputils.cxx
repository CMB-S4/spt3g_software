#include <pybindings.h>

#include <vector>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <G3Logging.h>
#include <G3Units.h>

#include <maps/maputils.h>
#include <maps/pointing.h>

#include <iostream>

#define COS cos
#define SIN sin
#define ATAN2 atan2

using namespace G3Units;

G3SkyMapPtr GetMaskMap(G3SkyMapConstPtr m, bool zero_nans, bool zero_infs)
{
	G3SkyMapPtr mask = m->Clone(false);

	for (size_t i = 0; i < m->size(); i++) {
		double v = m->at(i);
		if (zero_nans && v != v)
			continue;
		if (zero_infs && !std::isfinite(v))
			continue;
		if (v == 0)
			continue;
		(*mask)[i] = 1.0;
	}

	return mask;
}

void ApplyMask(G3SkyMapPtr m, G3SkyMapConstPtr mask, bool inverse)
{
	g3_assert(m->IsCompatible(*mask));

	for (size_t i = 0; i < m->size(); i++) {
		if (!(m->at(i)))
			continue;
		if (!!(mask->at(i)) == inverse)
			(*m)[i] = 0;
	}
}

void RemoveWeightsT(G3SkyMapPtr T, G3SkyMapWeightsConstPtr W, bool zero_nans)
{
	RemoveWeights(T, NULL, NULL, W, zero_nans);
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
				double c = m.Cond();
				empty = (c != c) || (c > 1e12);
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

void ApplyWeightsT(G3SkyMapPtr T, G3SkyMapWeightsConstPtr W)
{
	ApplyWeights(T, NULL, NULL, W);
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
	ra->units = G3Timestream::Angle;
	ra->pol_type = G3SkyMap::None;
	ra->SetPolConv(G3SkyMap::ConvNone);
	dec->weighted = false;
	dec->units = G3Timestream::Angle;
	dec->pol_type = G3SkyMap::None;
	dec->SetPolConv(G3SkyMap::ConvNone);

	return boost::python::make_tuple(ra, dec);
}


static double wrap_ra(double ra)
{
	static const double circ = 2 * M_PI * G3Units::rad;
	return fmod((ra < 0) ? (circ * (1. + ceilf(fabs(ra) / circ)) + ra) : ra, circ);
}


G3SkyMapPtr GetRaDecMask(G3SkyMapConstPtr m, double ra_left, double ra_right,
    double dec_bottom, double dec_top)
{
	G3SkyMapPtr mask = m->Clone(false);
	mask->weighted = false;
	mask->units = G3Timestream::None;
	mask->pol_type = G3SkyMap::None;
	mask->SetPolConv(G3SkyMap::ConvNone);

	ra_left = wrap_ra(ra_left);
	ra_right = wrap_ra(ra_right);

	for (size_t i = 0; i < m->size(); i++) {
		std::vector<double> radec = m->PixelToAngle(i);
		double ra = wrap_ra(radec[0]);
		if (ra_left < ra_right && (ra <= ra_left || ra >= ra_right))
			continue;
		if (ra_left >= ra_right && (ra <= ra_left && ra >= ra_right))
			continue;
		double dec = radec[1];
		if (dec <= dec_bottom || dec >= dec_top)
			continue;
		(*mask)[i] = 1.0;
	}

	return mask;
}


void FlattenPol(FlatSkyMapPtr Q, FlatSkyMapPtr U, G3SkyMapWeightsPtr W, double h, bool invert)
{
	if (U->GetPolConv() == G3SkyMap::ConvNone)
		log_warn("Missing pol_conv attribute for flatten_pol, assuming "
			 "U.pol_conv is set to IAU. This will raise an error "
			 "in the future.");

	if (!W)
		log_warn("Missing weights for flatten_pol");

	g3_assert(Q->IsCompatible(*U));
	g3_assert(Q->IsPolFlat() == U->IsPolFlat());
	FlatSkyMapPtr flatptr;
	if (!!W) {
		flatptr = boost::dynamic_pointer_cast<FlatSkyMap>(W->TQ);
		g3_assert(flatptr->IsCompatible(*Q));
		g3_assert(flatptr->IsPolFlat() == Q->IsPolFlat());
	}

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
		if (U->GetPolConv() == G3SkyMap::COSMO)
			rot *= -1.0;
		double cr = COS(rot);
		double sr = SIN(rot);

		(*Q)[i.first] = cr * q - sr * u;
		(*U)[i.first] = sr * q + cr * u;

		if (!W)
			continue;

		MuellerMatrix w = (*W)[i.first];

		double sr2 = 2 * sr * cr;
		double cr2 = 1 - 2 * sr * sr;
		double ws = (w.qq + w.uu) / 2.0;
		double wd = (w.qq - w.uu) / 2.0;
		double delta = wd * cr2 - w.qu * sr2;
		double tq = w.tq;
		double tu = w.tu;
		w.tq = cr * tq - sr * tu;
		w.tu = sr * tq + cr * tu;
		w.qq = ws + delta;
		w.uu = ws - delta;
		w.qu = wd * sr2 + w.qu * cr2;
	}

	Q->SetFlatPol(!invert);
	U->SetFlatPol(!invert);

	if (!!W) {
		flatptr = boost::dynamic_pointer_cast<FlatSkyMap>(W->TQ);
		flatptr->SetFlatPol(!invert);
		flatptr = boost::dynamic_pointer_cast<FlatSkyMap>(W->TU);
		flatptr->SetFlatPol(!invert);
		flatptr = boost::dynamic_pointer_cast<FlatSkyMap>(W->QQ);
		flatptr->SetFlatPol(!invert);
		flatptr = boost::dynamic_pointer_cast<FlatSkyMap>(W->QU);
		flatptr->SetFlatPol(!invert);
		flatptr = boost::dynamic_pointer_cast<FlatSkyMap>(W->UU);
		flatptr->SetFlatPol(!invert);
	}
}


void ReprojMap(G3SkyMapConstPtr in_map, G3SkyMapPtr out_map, int rebin, bool interp)
{
	bool rotate = false; // no transform
	quat q_rot; // quaternion for rotating from output to input coordinate system
	if (in_map->coord_ref != out_map->coord_ref &&
	    in_map->coord_ref != MapCoordReference::Local &&
	    out_map->coord_ref != MapCoordReference::Local) {
		rotate = true;
		q_rot = get_fk5_j2000_to_gal_quat();
		if (in_map->coord_ref == MapCoordReference::Equatorial)
			q_rot = 1. / q_rot;
	} else if (in_map->coord_ref != out_map->coord_ref) {
		log_fatal("Cannot convert input coord_ref %d to output coord_ref %d",
		    in_map->coord_ref, out_map->coord_ref);
	}

	if (out_map->pol_type != G3SkyMap::None && out_map->pol_type != in_map->pol_type) {
		log_fatal("Cannot convert input pol_type %d to output pol_type %d",
		    in_map->pol_type, out_map->pol_type);
	} else {
		out_map->pol_type = in_map->pol_type;
	}

	double s = 1.;
	if (out_map->GetPolConv() != G3SkyMap::ConvNone && in_map->GetPolConv() != G3SkyMap::ConvNone) {
		if (out_map->pol_type == G3SkyMap::U && out_map->GetPolConv() != in_map->GetPolConv())
			s = -1.;
	} else if (out_map->GetPolConv() == G3SkyMap::ConvNone) {
		out_map->SetPolConv(in_map->GetPolConv());
	}

	for (size_t i = 0; i < out_map->size(); i++) {
		double val = 0;
		if (rebin > 1) {
			std::vector<double> ra, dec;
			out_map->GetRebinAngles(i, rebin, ra, dec);
			for (size_t j = 0; j < ra.size(); j++) {
				if (rotate) {
					quat q = ang_to_quat(ra[j], dec[j]);
					q = q_rot * q / q_rot;
					quat_to_ang(q, ra[j], dec[j]);
				}
				if (interp)
					val += in_map->GetInterpValue(ra[j], dec[j]);
				else
					val += in_map->at(in_map->AngleToPixel(ra[j], dec[j]));
			}
			val /= ra.size();
		} else {
			std::vector<double> radec = out_map->PixelToAngle(i);
			if (rotate) {
				quat q = ang_to_quat(radec[0], radec[1]);
				q = q_rot * q / q_rot;
				quat_to_ang(q, radec[0], radec[1]);
			}
			if (interp)
				val = in_map->GetInterpValue(radec[0], radec[1]);
			else
				val = in_map->at(in_map->AngleToPixel(radec[0], radec[1]));
		}
		if (val != 0)
			(*out_map)[i] = s * val;
	}

	out_map->weighted = in_map->weighted;
	out_map->units = in_map->units;
}

// algorithm from https://www.johndcook.com/blog/skewness_kurtosis/
std::vector<double> GetMapStats(G3SkyMapConstPtr m, G3SkyMapConstPtr mask,
    int order, bool ignore_zeros, bool ignore_nans, bool ignore_infs)
{
	size_t n = 0;
	double m1 = 0;
	double m2 = 0;
	double m3 = 0;
	double m4 = 0;
	double a, b, c;

	for (size_t i = 0; i < m->size(); i++) {
		if (!!mask && !mask->at(i))
			continue;
		double v = m->at(i);
		if (ignore_zeros && v == 0)
			continue;
		if (ignore_nans && v != v)
			continue;
		if (ignore_infs && !std::isfinite(v))
			continue;

		n++;
		a = (v - m1) / n;
		if (order > 1) {
			b = a * a;
			c = b * n * (n - 1);
		}

		m1 += a;
		if (order > 3)
			m4 += c * b * (n * n - 3 * n + 3) + 6 * b * m2 - 4 * a * m3;
		if (order > 2)
			m3 += c * a * (n - 2) - 3 * a * m2;
		if (order > 1)
			m2 += c;
	}

	std::vector<double> out = {m1};

	if (order > 1)
		out.push_back(m2 / n);
	if (order > 2)
		out.push_back(sqrt((double)n) * m3/ pow(m2, 1.5));
	if (order > 3)
		out.push_back(n * m4 / (m2 * m2) - 3.0);

	return out;
}


double
GetMapMedian(G3SkyMapConstPtr m, G3SkyMapConstPtr mask, bool ignore_zeros, bool ignore_nans,
    bool ignore_infs)
{

	std::vector<double> data;
	size_t npix = mask ? mask->NpixNonZero() : (ignore_zeros ? m->NpixAllocated() : m->size());
	if (npix == 0)
		return 0;

	data.reserve(npix);

	for (size_t i = 0; i < m->size(); i++) {
		if (!!mask && !mask->at(i))
			continue;
		double v = m->at(i);
		if (ignore_zeros && v == 0)
			continue;
		if (ignore_nans && v != v)
			continue;
		if (ignore_infs && !std::isfinite(v))
			continue;
		data.push_back(v);
	}

	npix = data.size();
	size_t n = npix / 2;
	std::nth_element(data.begin(), data.begin() + n, data.end());

	// odd-length array
	if (npix % 2)
		return data[n];

	// even-length array
	double median = data[n];
	std::nth_element(data.begin(), data.begin() + n - 1, data.end());
	return (median + data[n - 1]) / 2.;
}


std::vector<double>
GetMapMinMax(G3SkyMapConstPtr m, G3SkyMapConstPtr mask, bool ignore_zeros, bool ignore_nans,
    bool ignore_infs)
{
	double min = 0.0 / 0.0;
	double max = 0.0 / 0.0;

	for (size_t i = 0; i < m->size(); i++) {
		if (!!mask && !mask->at(i))
			continue;
		double v = m->at(i);
		if (ignore_zeros && v == 0)
			continue;
		if (ignore_nans && v != v)
			continue;
		if (ignore_infs && !std::isfinite(v))
			continue;
		if (v > max || max != max)
			max = v;
		if (v < min || min != min)
			min = v;
	}

	return {min, max};
}


std::vector<double>
GetMapHist(G3SkyMapConstPtr m, const std::vector<double> &bin_edges, G3SkyMapConstPtr mask,
    bool ignore_zeros, bool ignore_nans, bool ignore_infs)
{

	g3_assert(std::is_sorted(bin_edges.begin(), bin_edges.end()));

	double bin_min = bin_edges[0];
	double bin_max = bin_edges[bin_edges.size() - 1];
	size_t nbins = bin_edges.size() - 1;
	double bin_width = (bin_max - bin_min) / (float) nbins;

	// determine whether bins are equally spaced
	bool equal_bins = true;
	for (size_t i = 1; i < bin_edges.size(); i++) {
		double bw = bin_edges[i] - bin_edges[i - 1];
		if (std::abs(bw - bin_width) > 1e-8) {
			equal_bins = false;
			break;
		}
	}

	std::vector<double> hist(nbins);

	for (size_t i = 0; i < m->size(); i++) {
		if (!!mask && !mask->at(i))
			continue;
		double v = m->at(i);
		if (ignore_zeros && v == 0)
			continue;
		if (ignore_nans && v != v)
			continue;
		if (ignore_infs && !std::isfinite(v))
			continue;
		if ((v < bin_min) || (v > bin_max))
			continue;

		size_t bin;
		if (v == bin_max) {
			// mimic numpy hist behavior, right-most bin edge is closed
			bin = nbins - 1;
		} else if (equal_bins) {
			// fast path for equally spaced bins
			bin = (size_t)floor((v - bin_min) / bin_width);
		} else {
			// log(N) search for bin number
			auto it = std::upper_bound(bin_edges.begin(), bin_edges.end(), v);
			bin = it - bin_edges.begin() - 1;
		}

		hist[bin] += 1;
	}

	return hist;
}


FlatSkyMapPtr
ConvolveMap(FlatSkyMapConstPtr map, FlatSkyMapConstPtr kernel)
{
	int xdim = map->shape()[0];
	int ydim = map->shape()[1];
	int nx = kernel->shape()[0];
	int ny = kernel->shape()[1];
	if ((nx % 2 == 0) || (ny % 2 == 0))
		log_fatal("Kernel must have odd map dimensions");

	FlatSkyMapPtr outmap = boost::dynamic_pointer_cast<FlatSkyMap>(map->Clone(false));
	if (map->IsDense())
		outmap->ConvertToDense();

	// loop over only non-zero kernel values
	std::vector<int> xk, yk;
	std::vector<double> vk;
	for (auto i: *kernel) {
		if (i.second == 0)
			continue;
		vk.push_back(i.second);
		xk.push_back(nx / 2 - (i.first % nx));
		yk.push_back(ny / 2 - (i.first / nx));
	}
	size_t nk = vk.size();

	for (size_t y = 0; y < ydim; y++) {
		for (size_t x = 0; x < xdim; x++) {
			double v = 0;
			for (size_t j = 0; j < nk; j++) {
				double m = map->at(x + xk[j], y + yk[j]);
				if (m != 0)
					v += vk[j] * m;
			}
			if (v != 0)
				(*outmap)(x, y) = v;
		}
	}

	return outmap;
}


static FlatSkyMapPtr
pyconvolve_map(FlatSkyMapConstPtr map, bp::object val)
{

	FlatSkyMapConstPtr kernel;
	if (bp::extract<FlatSkyMap>(val).check())
		kernel = bp::extract<FlatSkyMapConstPtr>(val)();
	else
		kernel = FlatSkyMapConstPtr(new FlatSkyMap(val, map->yres()));
	return ConvolveMap(map, kernel);
}


namespace bp = boost::python;
void maputils_pybindings(void){
	bp::def("get_mask_map", GetMaskMap,
		(bp::arg("map_in"), bp::arg("zero_nans")=false, bp::arg("zero_infs")=false),
		"Returns a map that is 1 where the input map is nonzero, and optionally "
		"non-nan and/or non-infinite");

	bp::def("apply_mask", ApplyMask,
		(bp::arg("map"), bp::arg("mask"), bp::arg("inverse")=false),
		"Apply a boolean mask to the input map in place, optionally inverting the mask values");

	bp::def("remove_weights_t", RemoveWeightsT,
		(bp::arg("T"), bp::arg("W"), bp::arg("zero_nans")=false),
		"Remove weights from unpolarized maps.	If zero_nans is true, empty pixels "
		"are skipped, and pixels with zero weight are set to 0 instead of nan.");

	bp::def("remove_weights", RemoveWeights,
		(bp::arg("T"), bp::arg("Q"), bp::arg("U"), bp::arg("W"), bp::arg("zero_nans")=false),
		"Remove weights from polarized maps.  If zero_nans is true, empty pixels "
		"are skipped, and pixels with zero weight are set to 0 instead of nan.");

	bp::def("apply_weights_t", ApplyWeightsT,
		(bp::arg("T"), bp::arg("W")),
		"Apply weights to unpolarized maps.");

	bp::def("apply_weights", ApplyWeights,
		(bp::arg("T"), bp::arg("Q"), bp::arg("U"), bp::arg("W")),
		"Apply weights to polarized maps.");

	bp::def("get_ra_dec_map", GetRaDecMap, (bp::arg("map_in")),
		"Returns maps of the ra and dec angles for each pixel in the input map");

	bp::def("get_ra_dec_mask", GetRaDecMask,
		(bp::arg("map_in"), bp::arg("ra_left"), bp::arg("ra_right"),
		 bp::arg("dec_bottom"), bp::arg("dec_top")),
		"Returns a mask that is nonzero for any pixels within the given ra and dec ranges");

	bp::def("flatten_pol", FlattenPol,
		(bp::arg("Q"), bp::arg("U"), bp::arg("W")=G3SkyMapWeightsPtr(),
		 bp::arg("h")=0.001, bp::arg("invert")=false),
		"For maps defined on the sphere the direction of the polarization angle is "
		"is defined relative to the direction of North.  When making maps we follow "
		"this definition.\n\nFor any flat sky estimators, the polarization angle is "
		"defined relative to the vertical axis.  For some map projections the "
		"direction of north is not the same as the vertical axis.  This function "
		"applies a rotation to the Q and U values to switch the curved sky Q/U "
		"definition to the flat sky Q/U definition.\n\nIf for whatever reason you "
		"want to reverse the process set the invert argument to True. Also applies "
		"the appropriate rotation to the Q and u elements of the associated weights.");

	bp::def("reproj_map", ReprojMap,
		(bp::arg("in_map"), bp::arg("out_map"), bp::arg("rebin")=1, bp::arg("interp")=false),
		"Reprojects the data from in_map onto out_map.  out_map can have a different "
		"projection, size, resolution, etc.  Optionally account for sub-pixel "
		"structure by setting rebin > 1 and/or enable bilinear interpolation of "
		"values from the input map by setting interp=True.  Use the maps' coord_ref "
		"attributes to rotate between Equatorial and Galactic coordinate systems.  "
		"Use the maps' pol_conv attributes to switch between COSMO and IAU "
		"polarization conventions.  If output attributes are not set, they will be "
		"copied from the input map.");

	bp::def("get_map_stats", GetMapStats,
		(bp::arg("map"), bp::arg("mask")=G3SkyMapConstPtr(), bp::arg("order")=2,
		 bp::arg("ignore_zeros")=false, bp::arg("ignore_nans")=false, bp::arg("ignore_infs")=false),
		"Computes moment statistics of the input map, optionally ignoring "
		"zero, nan and/or inf values in the map.  If order = 1, only the mean is "
		"returned.  If order = 2, 3 or 4 then the variance, skew and kurtosis "
		"are also included, respectively.  If a mask is supplied, then only "
		"the non-zero pixels in the mask are included.");

	bp::def("get_map_median", GetMapMedian,
		(bp::arg("map"), bp::arg("mask")=G3SkyMapConstPtr(),
		 bp::arg("ignore_zeros")=false, bp::arg("ignore_nans")=false, bp::arg("ignore_infs")=false),
		"Computes the median of the input map, optionally ignoring zero, nan and/or inf "
		"values in the map.  Requires making a copy of the data.  If a mask is "
		"supplied, then only the non-zero pixels in the mask are included.");

	bp::def("get_map_minmax", GetMapMinMax,
		(bp::arg("map"), bp::arg("mask")=G3SkyMapConstPtr(),
		 bp::arg("ignore_zeros")=false, bp::arg("ignore_nans")=false, bp::arg("ignore_infs")=false),
		"Computes the min and max of the input map, optionally ignoring "
		"zero, nan and/or inf values in the map.  If a mask is supplied, then "
		"only the non-zero pixels in the mask are included.");

	bp::def("get_map_hist", GetMapHist,
		(bp::arg("map"), bp::arg("bin_edges"), bp::arg("mask")=G3SkyMapConstPtr(),
		 bp::arg("ignore_zeros")=false, bp::arg("ignore_nans")=false, bp::arg("ignore_infs")=false),
		"Computes the histogram of the input map into bins defined by the array of "
		"bin edges, optionally ignoring zero, nan and/or inf values in the map.  "
		"If a mask is supplied, then only the non-zero pixels in the mask are included.");

	bp::def("convolve_map", pyconvolve_map, (bp::arg("map"), bp::arg("kernel")),
		"Convolve the input flat sky map with the given map-space kernel. The "
		"kernel must have odd dimensions and the same resolution as the map.");
}
