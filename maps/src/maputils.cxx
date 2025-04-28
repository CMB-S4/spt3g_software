#include <pybindings.h>

#include <vector>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <G3Logging.h>
#include <G3Units.h>

#include <maps/G3SkyMap.h>
#include <maps/FlatSkyMap.h>
#include <maps/pointing.h>

#include <iostream>

#define COS cos
#define SIN sin
#define ATAN2 atan2

using namespace G3Units;

void RemoveWeights(G3SkyMap &T, G3SkyMap &Q, G3SkyMap &U, const G3SkyMapWeights &W,
    bool zero_nans)
{
	g3_assert(W.IsPolarized());
	g3_assert(T.weighted);
	g3_assert(W.IsCongruent());
	g3_assert(T.IsCompatible(*(W.TT)));

	g3_assert(T.IsCompatible(Q));
	g3_assert(T.IsCompatible(U));
	g3_assert(Q.weighted);
	g3_assert(U.weighted);

	if (!zero_nans) {
		T.ConvertToDense();
		Q.ConvertToDense();
		U.ConvertToDense();
		for (size_t pix = 0; pix < T.size(); pix++) {
			StokesVector v(T[pix], Q[pix], U[pix]);
			v /= W.at(pix);
		}
	} else {
		for (size_t pix = 0; pix < W.size(); pix++) {
			double t = T.at(pix);
			MuellerMatrix m = W.at(pix);
			bool empty;

			// skip empty pixels
			double c = m.cond();
			empty = (c != c) || (c > 1e12);
			if (empty && t == 0 && Q.at(pix) == 0 && U.at(pix) == 0)
				continue;
			if (!empty)
				empty |= (m.det() == 0);

			// set bad pixels to 0
			if (empty) {
				T[pix] = 0;
				Q[pix] = 0;
				U[pix] = 0;
				continue;
			}

			// remove weights
			StokesVector v(T[pix], Q[pix], U[pix]);
			v /= m;
		}
	}

	T.weighted = false;
	Q.weighted = false;
	U.weighted = false;
}

void RemoveWeightsT(G3SkyMap &T, const G3SkyMapWeights &W, bool zero_nans)
{
	g3_assert(!W.IsPolarized());
	g3_assert(T.weighted);
	g3_assert(W.IsCongruent());
	g3_assert(T.IsCompatible(*(W.TT)));

	if (!zero_nans) {
		T.ConvertToDense();
		T /= *(W.TT);
	} else {
		for (size_t pix = 0; pix < W.size(); pix++) {
			double t = T.at(pix);
			MuellerMatrix m = W.at(pix);
			bool empty;

			// skip empty pixels
			empty = (m.tt == 0);
			if (empty && t == 0)
				continue;

			// set bad pixels to 0
			if (empty) {
				T[pix] = 0;
				continue;
			}

			// remove weights
			T[pix] /= (*(W.TT))[pix];
		}
	}

	T.weighted = false;
}

void ApplyWeights(G3SkyMap &T, G3SkyMap &Q, G3SkyMap &U, const G3SkyMapWeights &W)
{
	g3_assert(W.IsPolarized());
	g3_assert(!T.weighted);
	g3_assert(W.IsCongruent());
	g3_assert(T.IsCompatible(*(W.TT)));

	g3_assert(T.IsCompatible(Q));
	g3_assert(T.IsCompatible(U));
	g3_assert(!Q.weighted);
	g3_assert(!U.weighted);

	for (size_t pix = 0; pix < T.size(); pix++) {
		if (T.at(pix) == 0 && Q.at(pix) == 0 && U.at(pix) == 0)
			continue;
		StokesVector v(T[pix], Q[pix], U[pix]);
		v = W.at(pix) * v;
	}

	T.weighted = true;
	Q.weighted = true;
	U.weighted = true;
}

void ApplyWeightsT(G3SkyMap &T, const G3SkyMapWeights &W)
{
	g3_assert(!W.IsPolarized());
	g3_assert(!T.weighted);
	g3_assert(W.IsCongruent());
	g3_assert(T.IsCompatible(*(W.TT)));

	T *= *(W.TT);
	T.weighted = true;
}

py::tuple GetRaDecMap(const G3SkyMap &m)
{

	G3SkyMapPtr ra = m.Clone(false);
	G3SkyMapPtr dec = m.Clone(false);

	// These are going to be dense maps, so just start that way
	ra->ConvertToDense();
	dec->ConvertToDense();

	for (size_t i = 0; i < m.size(); i++) {
		std::vector<double> radec = m.PixelToAngle(i);
		(*ra)[i] = radec[0];
		(*dec)[i] = radec[1];
	}

	ra->weighted = false;
	ra->units = G3Timestream::Angle;
	ra->pol_type = G3SkyMap::None;
	ra->pol_conv = G3SkyMap::ConvNone;
	dec->weighted = false;
	dec->units = G3Timestream::Angle;
	dec->pol_type = G3SkyMap::None;
	dec->pol_conv = G3SkyMap::ConvNone;

	return py::make_tuple(ra, dec);
}


static double wrap_ra(double ra)
{
	static const double circ = 2 * M_PI * G3Units::rad;
	return fmod((ra < 0) ? (circ * (1. + ceilf(fabs(ra) / circ)) + ra) : ra, circ);
}


G3SkyMapMaskPtr GetRaDecMask(const G3SkyMap &m, double ra_left, double ra_right,
    double dec_bottom, double dec_top)
{
	G3SkyMapMaskPtr mask(new G3SkyMapMask(m));

	ra_left = wrap_ra(ra_left);
	ra_right = wrap_ra(ra_right);

	for (size_t i = 0; i < m.size(); i++) {
		std::vector<double> radec = m.PixelToAngle(i);
		double ra = wrap_ra(radec[0]);
		if (ra_left < ra_right && (ra <= ra_left || ra >= ra_right))
			continue;
		if (ra_left >= ra_right && (ra <= ra_left && ra >= ra_right))
			continue;
		double dec = radec[1];
		if (dec <= dec_bottom || dec >= dec_top)
			continue;
		(*mask)[i] = true;
	}

	return mask;
}


void FlattenPol(FlatSkyMap &Q, FlatSkyMap &U, G3SkyMapWeightsPtr W, double h, bool invert)
{
	if (!(U.IsPolarized()))
		log_warn("Missing pol_conv attribute for flatten_pol, assuming "
			 "U.pol_conv is set to IAU. This will raise an error "
			 "in the future.");

	g3_assert(Q.IsCompatible(U));
	g3_assert(Q.IsPolFlat() == U.IsPolFlat());
	FlatSkyMapPtr flatptr;
	if (!!W) {
		g3_assert(W->IsCompatible(Q));
		flatptr = std::dynamic_pointer_cast<FlatSkyMap>(W->TQ);
		g3_assert(flatptr->IsPolFlat() == Q.IsPolFlat());
	}

	if (Q.IsPolFlat() && !invert)
		return;
	if (!Q.IsPolFlat() && invert)
		return;

	for (auto i : Q) {
		double q = i.second;
		double u = U.at(i.first);
		if (q == 0 && u == 0)
			continue;

		std::vector<double> grad = Q.PixelToAngleGrad(i.first, h);
		double rot = ATAN2(-grad[0], grad[1]) + ATAN2(-grad[3], -grad[2]);
		if (invert)
			rot *= -1.0;
		if (U.pol_conv == G3SkyMap::COSMO)
			rot *= -1.0;
		double cr = COS(rot);
		double sr = SIN(rot);

		Q[i.first] = cr * q - sr * u;
		U[i.first] = sr * q + cr * u;

		if (!W)
			continue;

		MuellerMatrix w = (*W)[i.first];

		double sr2 = 2.0 * sr * cr;
		double cr2 = 1.0 - 2.0 * sr * sr;
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

	Q.SetFlatPol(!invert);
	U.SetFlatPol(!invert);

	if (!!W) {
		flatptr = std::dynamic_pointer_cast<FlatSkyMap>(W->TT);
		flatptr->SetFlatPol(!invert);
		flatptr = std::dynamic_pointer_cast<FlatSkyMap>(W->TQ);
		flatptr->SetFlatPol(!invert);
		flatptr = std::dynamic_pointer_cast<FlatSkyMap>(W->TU);
		flatptr->SetFlatPol(!invert);
		flatptr = std::dynamic_pointer_cast<FlatSkyMap>(W->QQ);
		flatptr->SetFlatPol(!invert);
		flatptr = std::dynamic_pointer_cast<FlatSkyMap>(W->QU);
		flatptr->SetFlatPol(!invert);
		flatptr = std::dynamic_pointer_cast<FlatSkyMap>(W->UU);
		flatptr->SetFlatPol(!invert);
	}
}


void ReprojMap(const G3SkyMap &in_map, G3SkyMap &out_map, int rebin, bool interp,
    G3SkyMapMaskConstPtr out_map_mask)
{
	bool rotate = false; // no transform
	Quat q_rot; // quaternion for rotating from output to input coordinate system
	if (in_map.coord_ref != out_map.coord_ref &&
	    in_map.coord_ref != MapCoordReference::Local &&
	    out_map.coord_ref != MapCoordReference::Local) {
		rotate = true;
		q_rot = get_fk5_j2000_to_gal_quat();
		if (in_map.coord_ref == MapCoordReference::Equatorial)
			q_rot = ~q_rot;
	} else if (in_map.coord_ref != out_map.coord_ref) {
		log_fatal("Cannot convert input coord_ref %d to output coord_ref %d",
		    in_map.coord_ref, out_map.coord_ref);
	}

	if (out_map.pol_type != G3SkyMap::None && out_map.pol_type != in_map.pol_type) {
		log_fatal("Cannot convert input pol_type %d to output pol_type %d",
		    in_map.pol_type, out_map.pol_type);
	} else {
		out_map.pol_type = in_map.pol_type;
	}

	double s = 1.;
	if (out_map.IsPolarized() && in_map.IsPolarized()) {
		if (out_map.pol_type == G3SkyMap::U && out_map.pol_conv != in_map.pol_conv)
			s = -1.;
	} else if (!(out_map.IsPolarized())) {
		out_map.pol_conv = in_map.pol_conv;
	}

	if (!!out_map_mask && !out_map_mask->IsCompatible(out_map))
		log_fatal("Mask is not compatible with output map");

	size_t stop = out_map.size();
	if (rebin > 1) {
		for (size_t i = 0; i < stop; i++) {
			if (!!out_map_mask && !out_map_mask->at(i)) {
				out_map[i] = 0;
				continue;
			}
			double val = 0;
			auto quats = out_map.GetRebinQuats(i, rebin);
			if (rotate)
				quats = q_rot * quats * ~q_rot;
			if (interp)
				for (size_t j = 0; j < quats.size(); j++)
					val += in_map.GetInterpValue(quats[j]);
			else
				for (size_t j = 0; j < quats.size(); j++)
					val += in_map.at(in_map.QuatToPixel(quats[j]));
			if (val != 0) {
				val /= quats.size();
				out_map[i] = s * val;
			}
		}
	} else {
		for (size_t i = 0; i < stop; i++) {
			if (!!out_map_mask && !out_map_mask->at(i)) {
				out_map[i] = 0;
				continue;
			}
			double val = 0;
			auto q = out_map.PixelToQuat(i);
			if (rotate)
				q = q_rot * q * ~q_rot;
			if (interp)
				val = in_map.GetInterpValue(q);
			else
				val = in_map.at(in_map.QuatToPixel(q));
			if (val != 0)
				out_map[i] = s * val;
		}
	}

	out_map.weighted = in_map.weighted;
	out_map.units = in_map.units;
}

// algorithm from https://www.johndcook.com/blog/skewness_kurtosis/
std::vector<double> GetMapMoments(const G3SkyMap &m, G3SkyMapMaskConstPtr mask,
    int order, bool ignore_zeros, bool ignore_nans, bool ignore_infs)
{
	size_t n = 0;
	double m1 = 0;
	double m2 = 0;
	double m3 = 0;
	double m4 = 0;
	double a, b, c;

	for (size_t i = 0; i < m.size(); i++) {
		if (!!mask && !mask->at(i))
			continue;
		double v = m.at(i);
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


std::vector<double>
GetMapHist(const G3SkyMap &m, const std::vector<double> &bin_edges, G3SkyMapMaskConstPtr mask,
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

	for (size_t i = 0; i < m.size(); i++) {
		if (!!mask && !mask->at(i))
			continue;
		double v = m.at(i);
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
ConvolveMap(const FlatSkyMap &map, const FlatSkyMap &kernel)
{
	size_t xdim = map.shape()[0];
	size_t ydim = map.shape()[1];
	size_t nx = kernel.shape()[0];
	size_t ny = kernel.shape()[1];
	if ((nx % 2 == 0) || (ny % 2 == 0))
		log_fatal("Kernel must have odd map dimensions");

	FlatSkyMapPtr outmap = std::dynamic_pointer_cast<FlatSkyMap>(map.Clone(false));
	if (map.IsDense())
		outmap->ConvertToDense();

	// loop over only non-zero kernel values
	std::vector<ssize_t> xk, yk;
	std::vector<double> vk;
	for (auto i: kernel) {
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
				double m = map.at(x + xk[j], y + yk[j]);
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
pyconvolve_map(const FlatSkyMap &map, const py::object &val)
{

	if (py::isinstance<FlatSkyMap>(val))
		return ConvolveMap(map, val.cast<const FlatSkyMap &>());

	// reach into python
	auto pykernel = py::type::of<FlatSkyMap>()(val, map.yres());
	return ConvolveMap(map, pykernel.cast<const FlatSkyMap &>());
}


G3SkyMapMaskPtr
MakePointSourceMask(const G3SkyMap &map, const std::vector<double> & ra,
    const std::vector<double> & dec, const std::vector<double> & radius)
{
	G3SkyMapMaskPtr mask(new G3SkyMapMask(map));
	g3_assert(ra.size() == dec.size());
	g3_assert(ra.size() == radius.size());

	for (size_t i = 0; i < ra.size(); i++) {
		auto pixels = map.QueryDisc(ra[i], dec[i], radius[i]);
		for (auto p: pixels)
			(*mask)[p] = true;
	}

	return mask;
}


PYBINDINGS("maps", scope)
{
	scope.def("remove_weights_t", RemoveWeightsT,
		py::arg("T"), py::arg("W"), py::arg("zero_nans")=false,
		"Remove weights from unpolarized maps.	If zero_nans is true, empty pixels "
		"are skipped, and pixels with zero weight are set to 0 instead of nan.");

	scope.def("remove_weights", RemoveWeights,
		py::arg("T"), py::arg("Q"), py::arg("U"), py::arg("W"), py::arg("zero_nans")=false,
		"Remove weights from polarized maps.  If zero_nans is true, empty pixels "
		"are skipped, and pixels with zero weight are set to 0 instead of nan.");

	scope.def("apply_weights_t", ApplyWeightsT,
		py::arg("T"), py::arg("W"),
		"Apply weights to unpolarized maps.");

	scope.def("apply_weights", ApplyWeights,
		py::arg("T"), py::arg("Q"), py::arg("U"), py::arg("W"),
		"Apply weights to polarized maps.");

	scope.def("get_ra_dec_map", GetRaDecMap, py::arg("map_in"),
		"Returns maps of the ra and dec angles for each pixel in the input map");

	scope.def("get_ra_dec_mask", GetRaDecMask,
		py::arg("map_in"), py::arg("ra_left"), py::arg("ra_right"),
		py::arg("dec_bottom"), py::arg("dec_top"),
		"Returns a mask that is nonzero for any pixels within the given ra and dec ranges");

	scope.def("flatten_pol", FlattenPol,
		py::arg("Q"), py::arg("U"), py::arg("W")=G3SkyMapWeightsPtr(),
		py::arg("h")=0.001, py::arg("invert")=false,
		"For maps defined on the sphere the direction of the polarization angle is "
		"is defined relative to the direction of North.  When making maps we follow "
		"this definition.\n\nFor any flat sky estimators, the polarization angle is "
		"defined relative to the vertical axis.  For some map projections the "
		"direction of north is not the same as the vertical axis.  This function "
		"applies a rotation to the Q and U values to switch the curved sky Q/U "
		"definition to the flat sky Q/U definition.\n\nIf for whatever reason you "
		"want to reverse the process set the invert argument to True. Also applies "
		"the appropriate rotation to the Q and u elements of the associated weights.");

	scope.def("reproj_map", ReprojMap,
		py::arg("in_map"), py::arg("out_map"), py::arg("rebin")=1, py::arg("interp")=false,
		py::arg("mask")=G3SkyMapMaskConstPtr(),
		"Reprojects the data from in_map onto out_map.  out_map can have a different "
		"projection, size, resolution, etc.  Optionally account for sub-pixel "
		"structure by setting rebin > 1 and/or enable bilinear interpolation of "
		"values from the input map by setting interp=True.  Use the maps' coord_ref "
		"attributes to rotate between Equatorial and Galactic coordinate systems.  "
		"Use the maps' pol_conv attributes to switch between COSMO and IAU "
		"polarization conventions.  If output attributes are not set, they will be "
		"copied from the input map. out_map_mask, if given, skip the unused pixels"
		"and set these pixels to 0.");

	scope.def("get_map_moments", GetMapMoments,
		py::arg("map"), py::arg("mask")=G3SkyMapMaskConstPtr(), py::arg("order")=2,
		py::arg("ignore_zeros")=false, py::arg("ignore_nans")=false, py::arg("ignore_infs")=false,
		"Computes moment statistics of the input map, optionally ignoring "
		"zero, nan and/or inf values in the map.  If order = 1, only the mean is "
		"returned.  If order = 2, 3 or 4 then the variance, skew and kurtosis "
		"are also included, respectively.  If a mask is supplied, then only "
		"the non-zero pixels in the mask are included.");

	scope.def("get_map_hist", GetMapHist,
		py::arg("map"), py::arg("bin_edges"), py::arg("mask")=G3SkyMapMaskConstPtr(),
		py::arg("ignore_zeros")=false, py::arg("ignore_nans")=false, py::arg("ignore_infs")=false,
		"Computes the histogram of the input map into bins defined by the array of "
		"bin edges, optionally ignoring zero, nan and/or inf values in the map.  "
		"If a mask is supplied, then only the non-zero pixels in the mask are included.");

	scope.def("convolve_map", pyconvolve_map, py::arg("map"), py::arg("kernel"),
		"Convolve the input flat sky map with the given map-space kernel. The "
		"kernel must have odd dimensions and the same resolution as the map.");

	scope.def("make_point_source_mask", MakePointSourceMask,
		py::arg("map"), py::arg("ra"), py::arg("dec"), py::arg("radius"),
		"Construct a mask from the input stub map with pixels within the given "
		"radius around each point source position set to 1.");
}
