#include <pybindings.h>
#include <serialization.h>
#include <coordinateutils/G3SkyMap.h>

namespace bp=boost::python;

template <class A> void
G3SkyMap::serialize(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("coord_ref", coord_ref);
	ar & cereal::make_nvp("units", units);
	if (v == 1) {
		std::vector<double> dat;
		ar & cereal::make_nvp("data", dat);
		
		uint32_t xpix, ypix;
		ar & cereal::make_nvp("xpix", xpix);
		ar & cereal::make_nvp("ypix", ypix);
		std::vector<size_t> dims;
		dims.push_back(xpix);
		dims.push_back(ypix);

		overflow = dat[dat.size() - 1];
		dat.resize(dat.size() - 1);
		init_from_v1_data(dims, dat);
	} else {
		ar & cereal::make_nvp("overflow", overflow);
	}

	ar & cereal::make_nvp("pol_type", pol_type);
	ar & cereal::make_nvp("is_weighted", is_weighted);
}

template <class A> void
G3SkyMapWeights::serialize(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("TT", TT);
	ar & cereal::make_nvp("TQ", TQ);
	ar & cereal::make_nvp("TU", TU);
	ar & cereal::make_nvp("QQ", QQ);
	ar & cereal::make_nvp("QU", QU);
	ar & cereal::make_nvp("UU", UU);
	if (v > 1) {
		ar & cereal::make_nvp("weight_type", weight_type);
	}
}

template <class A> void
G3SkyMapWithWeights::serialize(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("map_id", map_id);
	ar & cereal::make_nvp("T", T);
        ar & cereal::make_nvp("Q", Q);
        ar & cereal::make_nvp("U", U);
        ar & cereal::make_nvp("weights", weights);
}

G3_SERIALIZABLE_CODE(G3SkyMap);
G3_SERIALIZABLE_CODE(G3SkyMapWeights);

G3SkyMapWeights::G3SkyMapWeights(G3SkyMapConstPtr ref, 
  G3SkyMapWeights::WeightType wt = G3SkyMapWeights::Wpol) :
   TT(ref->Clone(false)), TQ(ref->Clone(false)), TU(ref->Clone(false)),
   QQ(ref->Clone(false)), QU(ref->Clone(false)), UU(ref->Clone(false)),
   weight_type(wt)
{
}

G3SkyMapWeights::G3SkyMapWeights(const G3SkyMapWeights &r) :
    TT(r.TT->Clone(true)), TQ(r.TQ->Clone(true)), TU(r.TU->Clone(true)),
    QQ(r.QQ->Clone(true)), QU(r.QU->Clone(true)), UU(r.UU->Clone(true)),
    weight_type(r.weight_type)
{
}

std::vector<int> G3SkyMap::angles_to_pixels(const std::vector<double> & alphas,
    const std::vector<double> & deltas) const
{
	std::vector<int> pixels(alphas.size());

	for (size_t i = 0; i < alphas.size(); i++) {
		pixels[i] = angle_to_pixel(alphas[i], deltas[i]);
	}

	return pixels;
}

void G3SkyMap::pixels_to_angles(const std::vector<int> & pixels,
    std::vector<double> & alphas, std::vector<double> & deltas) const
{
	if (alphas.size() != pixels.size()) {
		alphas = std::vector<double>(pixels.size());
	}
	if (deltas.size() != pixels.size()) {
		deltas = std::vector<double>(pixels.size());
	}

	for (size_t i = 0; i < pixels.size(); i++) {
		std::vector<double> ang;
		ang = pixel_to_angle(pixels[i]);
		alphas[i] = ang[0];
		deltas[i] = ang[1];
	}
}

size_t G3SkyMap::size() const
{
	size_t s = 1;
	for (size_t i : shape())
		s *= i;
	return s;
}

G3SkyMap &G3SkyMap::operator+=(const G3SkyMap & rhs)
{
	assert(IsCompatible(rhs));
	for (size_t i = 0; i < rhs.size(); i++)
		(*this)[i] += rhs[i];
	return *this;
}

G3SkyMap &G3SkyMap::operator+=(double rhs)
{
	for (size_t i = 0; i < size(); i++)
		(*this)[i] += rhs;
	return *this;
}

G3SkyMap &G3SkyMap::operator-=(const G3SkyMap &rhs)
{
	assert(IsCompatible(rhs));
	for (size_t i = 0; i < rhs.size(); i++)
		(*this)[i] -= rhs[i];
	return *this;
}

G3SkyMap &G3SkyMap::operator-=(double rhs)
{
	for (size_t i = 0; i < size(); i++)
		(*this)[i] -= rhs;
	return *this;
}

G3SkyMap &G3SkyMap::operator*=(const G3SkyMap &rhs)
{
	assert(IsCompatible(rhs));
	for (size_t i = 0; i < rhs.size(); i++)
		(*this)[i] *= rhs[i];
	return *this;
}

G3SkyMap &G3SkyMap::operator*=(double rhs)
{
	for (size_t i = 0; i < size(); i++)
		(*this)[i] *= rhs;
	return *this;
}

G3SkyMap &G3SkyMap::operator/=(const G3SkyMap &rhs)
{
	assert(IsCompatible(rhs));
	for (size_t i = 0; i < rhs.size(); i++)
		(*this)[i] /= rhs[i];
	return *this;
}

G3SkyMap &G3SkyMap::operator/=(double rhs)
{
	for (size_t i = 0; i < size(); i++)
		(*this)[i] /= rhs;
	return *this;
}

double G3SkyMap::get_interp_precalc(const std::vector<long> & pix,
    const std::vector<double> & weight) const
{
	double outval = 0;
	for (size_t i = 0; i < pix.size(); i++) {
		outval += (*this)[pix[i]] * weight[i];
	}
	return outval;
}

double G3SkyMap::get_interp_value(double alpha, double delta) const
{
	std::vector<long> pix;
	std::vector<double> weight;
	get_interp_pixels_weights(alpha, delta, pix, weight);
	return get_interp_precalc(pix, weight);
}

std::vector<double> G3SkyMap::get_interp_values(const std::vector<double> & alphas,
    const std::vector<double> & deltas) const
{
	std::vector<double> outvals(alphas.size());

	for (size_t i = 0; i < alphas.size(); i++) {
		outvals[i] = get_interp_value(alphas[i], deltas[i]);
	}

	return outvals;
}

static double
skymap_getitem(const G3SkyMap &skymap, int i)
{

	if (i < 0)
		i = skymap.size() + i;
	if (size_t(i) >= skymap.size()) {
		PyErr_SetString(PyExc_IndexError, "Index out of range");
		bp::throw_error_already_set();
	}

	return skymap[i];
}

static void
skymap_setitem(G3SkyMap &skymap, int i, double val)
{

	if (i < 0)
		i = skymap.size() + i;
	if (size_t(i) >= skymap.size()) {
		PyErr_SetString(PyExc_IndexError, "Index out of range");
		bp::throw_error_already_set();
	}

	skymap[i] = val;
}

static bp::tuple
skymap_shape(G3SkyMap &skymap)
{
	// Swap to match numpy's convention for shape()
	std::vector<size_t> shape = skymap.shape();
	std::vector<uint64_t> pyshape;
	for (ssize_t i = shape.size() - 1; i >= 0; i--)
		pyshape.push_back(shape[i]);
	
	return bp::tuple(pyshape);
}

static G3SkyMapPtr
skymap_copy(G3SkyMap &r)
{
	return r.Clone(true);
}


StokesVector & StokesVector::operator /=(const MuellerMatrix &r)
{
	double det = r.det();
	if (r.tt == 0 || det < 1e-12) {
		if (det < 1e-12 && r.tt != 0) {
			log_trace("Singular matrix found when inverting!  Det is %lE\n", det);
		}
		t = 0;
		q = 0;
		u = 0;
		return *this;
	}

	double t_ = (t * (r.qq * r.uu - r.qu * r.qu) +
		     -q * (r.tq * r.uu - r.tu * r.qu) +
		     u * (r.tq * r.qu - r.tu * r.qq)) / det;
	double q_ = (-t * (r.tq * r.uu - r.qu * r.tu) +
		     q * (r.tt * r.uu - r.tu * r.tu) +
		     -u * (r.tt * r.qu - r.tq * r.tu)) / det;
	double u_ = (t * (r.tq * r.qu - r.tu * r.qq) +
		     -q * (r.tt * r.qu - r.tu * r.tq) +
		     u * (r.tt * r.qq - r.tq * r.tq)) / det;
	t = t_;
	q = q_;
	u = u_;
	return *this;
}

StokesVector StokesVector::operator /(const MuellerMatrix &r) const
{
	StokesVector v;
	v += *this;
	v /= r;
	return v;
}

void G3SkyMapWithWeights::ApplyWeights(G3SkyMapWeightsPtr w)
{
	g3_assert(!IsWeighted());
	g3_assert(T->IsCompatible(*(w->TT)));

	for (size_t pix = 0; pix < T->size(); pix++) {
		// XXX implicitly densifies maps
		(*this)[pix] = (*w)[pix] * (*this)[pix];
	}

	// Store pointer to weights here
	weights = w;
}

void G3SkyMapWithWeights::RemoveWeights()
{
	g3_assert(IsWeighted());

	for (size_t pix = 0; pix < T->size(); pix++) {
		// XXX implicitly densifies maps
		(*this)[pix] /= (*weights)[pix];
	}

	// Remove pointer to weights
	weights.reset();
}

PYBINDINGS("coordinateutils") {
	bp::enum_<MapCoordReference>("MapCoordReference")
	    .value("Local", Local)
	    .value("Equatorial", Equatorial)
	    .value("Galactic", Galactic)
	;

	bp::enum_<G3SkyMap::MapPolType>("MapPolType")
	    .value("T", G3SkyMap::T)
	    .value("Q", G3SkyMap::Q)
	    .value("U", G3SkyMap::U)
	    .value("none", G3SkyMap::None) // "None" is reserved in python
	;
	enum_none_converter::from_python<G3SkyMap::MapPolType>();

	bp::enum_<G3SkyMapWeights::WeightType>("WeightType")
	    .value("Wpol", G3SkyMapWeights::Wpol)
	    .value("Wunpol", G3SkyMapWeights::Wunpol)
	    .value("None", G3SkyMapWeights::None)
	;
	enum_none_converter::from_python<G3SkyMapWeights::WeightType>();

	bp::class_<G3SkyMap, boost::noncopyable,
	  G3SkyMapPtr>("G3SkyMap",
	  "Base class for 1- and 2-D skymaps of various projections. Usually "
	  "you want a subclass of this (e.g. FlatSkyMap) rather than using it "
	  "directly.", bp::no_init)
	    .def_readwrite("coord_ref", &G3SkyMap::coord_ref,
	      "Coordinate system (coordinateutils.MapCoordReference) of the map (e.g. "
	      "Galactic, Equatorial, etc.)")
	    .def_readwrite("pol_type", &G3SkyMap::pol_type,
	      "Polarization type (coordinateutils.MapPolType) of the map "
	      "(e.g. coordinateutils.MapPolType.Q)")
	    .def_readwrite("units", &G3SkyMap::units,
	      "Unit class (core.G3TimestreamUnits) of the map (e.g. "
	      "core.G3TimestreamUnits.Tcmb). Within each unit class, further "
	      "conversions, for example from K to uK, should use core.G3Units.")
	    .def_readwrite("is_weighted", &G3SkyMap::is_weighted,
	      "True if map is multiplied by weights")
	    .add_property("size", &G3SkyMap::size, "Number of pixels in map")
	    .def("__len__", &G3SkyMap::size, "Number of pixels in map")
	    .add_property("shape", &skymap_shape, "Shape of map")
	    .add_property("npix_allocated", &G3SkyMap::npix_allocated,
	      "Number of pixels in map currently stored in memory")
	    .def_readwrite("overflow", &G3SkyMap::overflow,
              "Combined value of data processed by "
	      "the map maker but outside of the map area")
	    .def("__getitem__", &skymap_getitem)
	    .def("__setitem__", &skymap_setitem)
	    .def("__copy__", &skymap_copy)
	    .def("__deepcopy__", &skymap_copy)
	    .def("Clone", &G3SkyMap::Clone,
	      ((bp::arg("copy_data")=true),
	       "Return a map of the same type, populated with a copy of the data "
	       "if the argument is true (default), empty otherwise."))
	    .def("IsCompatible", &G3SkyMap::IsCompatible,
	      "Returns true if the input argument is a map with matching dimensions "
	      "and boundaries on the sky.")

	    .def("angles_to_pixels", &G3SkyMap::angles_to_pixels,
	      (bp::arg("alphas"), bp::arg("deltas")),
	       "Compute the 1D pixel location for each of the sky coordinates")
	    .def("pixels_to_angles", &G3SkyMap::pixels_to_angles,
	      (bp::arg("pixels"), bp::arg("alphas"), bp::arg("deltas")),
	       "Compute the sky coordinates of each of the given 1D pixels")
	    .def("pixel_to_angle", 
	      (std::vector<double> (G3SkyMap::*)(size_t) const) 
	      &G3SkyMap::pixel_to_angle, bp::arg("pixel"),
	      "Compute the sky coordinates of the given 1D pixel")
	    .def("angle_to_pixel", &G3SkyMap::angle_to_pixel,
	      (bp::arg("alpha"), bp::arg("delta")),
	       "Compute the 1D pixel location of the given sky coordinates.")

	    .def("get_interp_values", &G3SkyMap::get_interp_values,
	      (bp::arg("alphas"), bp::arg("deltas")),
	       "Return the values at each of the input coordinate locations. "
	       "Computes each value using bilinear interpolation over the "
	       "map pixels.")

	    .def("rebin", &G3SkyMap::rebin, bp::arg("scale"),
	      "Rebin the map into larger pixels by averaging scale-x-scale "
	      "blocks of pixels together.  Returns a new map object. "
	      "Map dimensions must be a multiple of the rebinning scale.")

	    .def(bp::self += bp::self)
	    .def(bp::self *= bp::self)
	    .def(bp::self -= bp::self)
	    .def(bp::self /= bp::self)
	    .def(bp::self += double())
	    .def(bp::self *= double())
	    .def(bp::self -= double())
	    .def(bp::self /= double())
	;

	EXPORT_FRAMEOBJECT(G3SkyMapWeights, init<>(), "generic sky weight")
	    .def(bp::init<G3SkyMapConstPtr, G3SkyMapWeights::WeightType>(
	      (bp::arg("skymap"),
               bp::arg("weight_type") = G3SkyMapWeights::Wpol)))
	    .def_readwrite("TT",&G3SkyMapWeights::TT)
	    .def_readwrite("TQ",&G3SkyMapWeights::TQ)
	    .def_readwrite("TU",&G3SkyMapWeights::TU)
	    .def_readwrite("QQ",&G3SkyMapWeights::QQ)
	    .def_readwrite("QU",&G3SkyMapWeights::QU)
	    .def_readwrite("UU",&G3SkyMapWeights::UU)
	    .def_readwrite("weight_type", &G3SkyMapWeights::weight_type)
	    
	    .def("Clone", &G3SkyMapWeights::Clone)
	;
	register_pointer_conversions<G3SkyMapWeights>();

#if 0
	EXPORT_FRAMEOBJECT(G3SkyMapWithWeights, init<>(), "Container for (potentially) polarized maps and weights")
	    .def(bp::init<G3SkyMapPtr, G3SkyMapWeightsPtr>(
	      (bp::arg("Tmap"),
               bp::arg("Tweights") = G3SkyMapWeightsPtr())))
	    .def(bp::init<G3SkyMapPtr, G3SkyMapPtr, G3SkyMapPtr, G3SkyMapWeightsPtr>(
	      (bp::arg("T"), bp::arg("Q"), bp::arg("U"),
               bp::arg("Tweights") = G3SkyMapWeightsPtr())))
	    .def_readwrite("T",&G3SkyMapWithWeights::T)
	    .def_readwrite("Q",&G3SkyMapWithWeights::Q)
	    .def_readwrite("U",&G3SkyMapWithWeights::U)
	    .def_readwrite("weights",&G3SkyMapWithWeights::weights)
	    .def_readwrite("map_id",&G3SkyMapWithWeights::map_id)
	    .add_property("weighted",&G3SkyMapWithWeights::IsWeighted)
	    .add_property("polarized",&G3SkyMapWithWeights::IsPolarized)
	    .def("remove_weights",&G3SkyMapWithWeights::RemoveWeights)
	    .def("apply_weights",&G3SkyMapWithWeights::ApplyWeights)
	    //.def("rebin", &G3SkyMapWithWeights::rebin) // XXX
	    
	    .def("Clone", &G3SkyMapWithWeights::Clone)
	;
	register_pointer_conversions<G3SkyMapWithWeights>();
#endif

}

