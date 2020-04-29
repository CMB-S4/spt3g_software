#include <pybindings.h>
#include <serialization.h>
#include <G3Units.h>
#include <maps/G3SkyMap.h>
#include <maps/FlatSkyMap.h>
#include <maps/HealpixSkyMap.h>

namespace bp=boost::python;

#define ACOS acos
#define COS cos

template <class A> void
G3SkyMap::serialize(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("coord_ref", coord_ref);
	ar & cereal::make_nvp("units", units);
	if (v == 1) {
		// Old versions had the buffer here instead of
		// in the derived class.
		std::vector<double> dat;
		ar & cereal::make_nvp("data", dat);
		
		uint32_t xpix, ypix;
		ar & cereal::make_nvp("xpix", xpix);
		ar & cereal::make_nvp("ypix", ypix);
		std::vector<size_t> dims;
		dims.push_back(xpix);
		dims.push_back(ypix);

		if (dat.size() > 0) {
			overflow = dat[dat.size() - 1];
			dat.resize(dat.size() - 1);
		} else {
			overflow = 0;
		}
		InitFromV1Data(dims, dat);
	} else {
		ar & cereal::make_nvp("overflow", overflow);
	}

	ar & cereal::make_nvp("pol_type", pol_type);
	ar & cereal::make_nvp("weighted", weighted);

	if (v > 2) {
		ar & cereal::make_nvp("pol_conv", pol_conv_);
	} else {
		pol_conv_ = ConvNone;
	}

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

	if (v == 2) {
		enum WeightType {
			Wpol = 3,
			Wunpol = 4,
			None = 5
		};
		WeightType weight_type = None;
		ar & cereal::make_nvp("weight_type", weight_type);
		if (weight_type == Wunpol) {
			TQ.reset();
			TU.reset();
			QQ.reset();
			QU.reset();
			UU.reset();
		}
	}
}

G3_SERIALIZABLE_CODE(G3SkyMap);
G3_SERIALIZABLE_CODE(G3SkyMapWeights);

G3SkyMap::G3SkyMap(MapCoordReference coords, bool weighted,
    G3Timestream::TimestreamUnits units, MapPolType pol_type,
    MapPolConv pol_conv) :
    coord_ref(coords), units(units), pol_type(pol_type),
    weighted(weighted), overflow(0), pol_conv_(pol_conv)
{
	if (pol_type == U && pol_conv == ConvNone)
		log_warn("Map object has pol_type U and unknown pol_conv. "
			 "Set the pol_conv attribute to IAU or COSMO.");
}

G3SkyMapWeights::G3SkyMapWeights(G3SkyMapConstPtr ref, bool polarized) :
    TT(ref->Clone(false)), TQ(polarized ? ref->Clone(false) : NULL),
    TU(polarized ? ref->Clone(false) : NULL),
    QQ(polarized ? ref->Clone(false) : NULL),
    QU(polarized ? ref->Clone(false) : NULL),
    UU(polarized ? ref->Clone(false) : NULL)
{
}

G3SkyMapWeights::G3SkyMapWeights(const G3SkyMapWeights &r) :
    TT(r.TT->Clone(true)), TQ(!r.TQ ? NULL : r.TQ->Clone(true)),
    TU(!r.TU ? NULL : r.TU->Clone(true)),
    QQ(!r.QQ ? NULL : r.QQ->Clone(true)),
    QU(!r.QU ? NULL : r.QU->Clone(true)),
    UU(!r.UU ? NULL : r.UU->Clone(true))
{
}

std::vector<int> G3SkyMap::AnglesToPixels(const std::vector<double> & alphas,
    const std::vector<double> & deltas) const
{
	std::vector<int> pixels(alphas.size());

	for (size_t i = 0; i < alphas.size(); i++) {
		pixels[i] = AngleToPixel(alphas[i], deltas[i]);
	}

	return pixels;
}

void G3SkyMap::PixelsToAngles(const std::vector<int> & pixels,
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
		ang = PixelToAngle(pixels[i]);
		alphas[i] = ang[0];
		deltas[i] = ang[1];
	}
}

static boost::python::tuple
skymap_pixels_to_angles(const G3SkyMap & skymap, const std::vector<int> & pixels)
{
	std::vector<double> alphas, deltas;
	skymap.PixelsToAngles(pixels, alphas, deltas);

	return boost::python::make_tuple(alphas, deltas);
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
	g3_assert(IsCompatible(rhs));
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
	g3_assert(IsCompatible(rhs));
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
	g3_assert(IsCompatible(rhs));
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
	g3_assert(IsCompatible(rhs));
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

double G3SkyMap::GetInterpPrecalc(const std::vector<long> & pix,
    const std::vector<double> & weight) const
{
	double outval = 0;
	for (size_t i = 0; i < pix.size(); i++) {
		outval += this->at(pix[i]) * weight[i];
	}
	return outval;
}

double G3SkyMap::GetInterpValue(double alpha, double delta) const
{
	std::vector<long> pix;
	std::vector<double> weight;
	GetInterpPixelsWeights(alpha, delta, pix, weight);
	return GetInterpPrecalc(pix, weight);
}

std::vector<double> G3SkyMap::GetInterpValues(const std::vector<double> & alphas,
    const std::vector<double> & deltas) const
{
	std::vector<double> outvals(alphas.size());

	for (size_t i = 0; i < alphas.size(); i++) {
		outvals[i] = GetInterpValue(alphas[i], deltas[i]);
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

	return skymap.at(i);
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

void
G3SkyMap::SetPolConv(G3SkyMap::MapPolConv pol_conv)
{
	if (pol_type == G3SkyMap::U && pol_conv == G3SkyMap::ConvNone)
		log_warn("Map object has pol_type U and unknown pol_conv. "
			 "Set the pol_conv attribute to IAU or COSMO.");

	if (pol_type != G3SkyMap::U ||
	    pol_conv == G3SkyMap::ConvNone ||
	    pol_conv_ == G3SkyMap::ConvNone) {
		pol_conv_ = pol_conv;
		return;
	}

	if (pol_conv != pol_conv_)
		(*this) *= -1;

	pol_conv_ = pol_conv;
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

static bp::tuple
skymapweights_shape(G3SkyMapWeights &weights)
{
	return skymap_shape(*(weights.TT));
}

static G3SkyMapPtr
skymap_copy(G3SkyMap &r)
{
	return r.Clone(true);
}

static G3SkyMapWeightsPtr
skymapweights_copy(G3SkyMapWeights &r)
{
	return r.Clone(true);
}

#define skymap_pynoninplace(name, oper, rhs_type) \
static G3SkyMapPtr \
pyskymap_##name(const G3SkyMap &a, const rhs_type b) \
{ \
	G3SkyMapPtr rv = a.Clone(true); \
	(*rv) oper b; \
	return rv; \
}

skymap_pynoninplace(add, +=, G3SkyMap &)
skymap_pynoninplace(addd, +=, double)
skymap_pynoninplace(sub, -=, G3SkyMap &)
skymap_pynoninplace(subd, -=, double)
skymap_pynoninplace(mult, *=, G3SkyMap &)
skymap_pynoninplace(multd, *=, double)
skymap_pynoninplace(div, /=, G3SkyMap &)
skymap_pynoninplace(divd, /=, double)

static G3SkyMapPtr
pyskymap_rsubd(const G3SkyMap &a, const double b)
{
	G3SkyMapPtr rv = pyskymap_subd(a, b);
	(*rv) *= -1;
	return rv;
}

static G3SkyMapPtr
pyskymap_rdivd(const G3SkyMap &a, const double b)
{
	// XXX this is pretty inefficient
	G3SkyMapPtr rv = a.Clone(false);
	(*rv) += b;
	(*rv) /= a;
	return rv;
}

static G3SkyMapPtr
pyskymap_neg(G3SkyMap &a)
{
	G3SkyMapPtr rv = a.Clone(true);
	(*rv) *= -1;
	return rv;
}

G3SkyMapWeights &
G3SkyMapWeights::operator+=(const G3SkyMapWeights &rhs)
{
	g3_assert(IsPolarized() == rhs.IsPolarized());

	*TT += *(rhs.TT);

	if (IsPolarized()) {
		*TQ += *rhs.TQ;
		*TU += *rhs.TU;
		*QQ += *rhs.QQ;
		*QU += *rhs.QU;
		*UU += *rhs.UU;
	}

	return *this;
}

#define skymapweights_inplace(op, rhs_type) \
G3SkyMapWeights &G3SkyMapWeights::operator op(rhs_type rhs) \
{ \
	*TT op rhs; \
	if (IsPolarized()) { \
		*TQ op rhs; \
		*TU op rhs; \
		*QQ op rhs; \
		*QU op rhs; \
		*UU op rhs; \
	} \
	return *this; \
}

skymapweights_inplace(*=, const G3SkyMap &);
skymapweights_inplace(*=, double);
skymapweights_inplace(/=, double);

#define skymapweights_pynoninplace(name, oper, rhs_type) \
static G3SkyMapWeightsPtr \
pyskymapweights_##name(const G3SkyMapWeights &a, const rhs_type b) \
{ \
	G3SkyMapWeightsPtr rv = a.Clone(true); \
	(*rv) oper b; \
	return rv; \
}

skymapweights_pynoninplace(add, +=, G3SkyMapWeights &);
skymapweights_pynoninplace(multm, *=, G3SkyMap &);
skymapweights_pynoninplace(multd, *=, double);
skymapweights_pynoninplace(divd, /=, double);

MuellerMatrix MuellerMatrix::Inv() const
{
	MuellerMatrix m;
	double d = Det();
	if (tt == 0 || d < 1e-12) {
		if (d < 1e-12 && tt != 0)
			log_trace("Singular matrix found when inverting!  Det is %lE\n", d);
		m.tt = m.tq = m.tu = m.qq = m.qu = m.uu = 0.0 / 0.0;
		return m;
	}

	m.tt = (qq * uu - qu * qu) / d;
	m.tq = (tu * qu - tq * uu) / d;
	m.tu = (tq * qu - tu * qq) / d;
	m.qq = (tt * uu - tu * tu) / d;
	m.qu = (tq * tu - tt * qu) / d;
	m.uu = (tt * qq - tq * tq) / d;

	return m;
}

double MuellerMatrix::Cond() const
{
	// Compute eigenvalues of a symmetrix 3x3 matrix
	// See https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3%C3%973_matrices
	double lmax, lmin;

	double p1 = tq * tq + tu * tu + qu * qu;
	if (p1 < 1e-12) {
		// matrix is empty
		if ((tt + qq + uu) == 0)
			return 0.0 / 0.0;

		// matrix is diagonal
		lmax = lmin = tt;
		if (qq > lmax)
			lmax = qq;
		if (qq < lmin)
			lmin = qq;
		if (uu > lmax)
			lmax = uu;
		if (uu < lmin)
			lmin = uu;
		return lmax / lmin;
	}

	double q = (tt + qq + uu) / 3.;
	double ttq = tt - q;
	double qqq = qq - q;
	double uuq = uu - q;
	double p = ttq * ttq + qqq * qqq + uuq * uuq + 2. * p1;
	p = sqrt(p / 6.0);

	MuellerMatrix B = *this;
	B.tt = ttq;
	B.qq = qqq;
	B.uu = uuq;
	B /= p;
	double r = B.Det() / 2.0;

	double phi;
	if (r <= -1)
		phi = 60 * G3Units::deg;
	else if (r >= 1)
		phi = 0;
	else
		phi = ACOS(r) / 3.0;

	lmax = q + 2. * p * COS(phi);
	lmin = q + 2. * p * COS(phi + 120 * G3Units::deg);
	return lmax / lmin;
}

StokesVector & StokesVector::operator /=(const MuellerMatrix &r)
{
	MuellerMatrix ir = r.Inv();
	if (ir.tt != ir.tt)
		t = q = u = 0.0 / 0.0;
	else
		(*this) = ir * (*this);
	return *this;
}

StokesVector StokesVector::operator /(const MuellerMatrix &r) const
{
	StokesVector v;
	v += *this;
	v /= r;
	return v;
}

G3SkyMapPtr G3SkyMapWeights::Det() const
{
	G3SkyMapPtr D = TT->Clone(false);

	for (size_t i = 0; i < TT->size(); i++) {
		double det = this->at(i).Det();
		if (det != 0)
			(*D)[i] = det;
	}

	return D;
}

G3SkyMapPtr G3SkyMapWeights::Cond() const
{
	G3SkyMapPtr C = TT->Clone(false);
	C->ConvertToDense();

	for (size_t i = 0; i < TT->size(); i++)
		(*C)[i] = this->at(i).Cond();

	return C;
}

G3SkyMapWeightsPtr G3SkyMapWeights::Rebin(size_t scale) const
{
	g3_assert(IsCongruent());
	G3SkyMapWeightsPtr out(new G3SkyMapWeights());

	out->TT = TT->Rebin(scale, false);
	if (IsPolarized()) {
		out->TQ = TQ->Rebin(scale, false);
		out->TU = TU->Rebin(scale, false);
		out->QQ = QQ->Rebin(scale, false);
		out->QU = QU->Rebin(scale, false);
		out->UU = UU->Rebin(scale, false);
	} else {
		out->TQ = NULL;
		out->TU = NULL;
		out->QQ = NULL;
		out->QU = NULL;
		out->UU = NULL;
	}

	return out;
}

void G3SkyMapWeights::Compact(bool zero_nans)
{
	g3_assert(IsCongruent());

	TT->Compact(zero_nans);
	if (IsPolarized()) {
		TQ->Compact(zero_nans);
		TU->Compact(zero_nans);
		QQ->Compact(zero_nans);
		QU->Compact(zero_nans);
		UU->Compact(zero_nans);
	}
}

PYBINDINGS("maps") {
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

	bp::enum_<G3SkyMap::MapPolConv>("MapPolConv")
	    .value("IAU", G3SkyMap::IAU)
	    .value("COSMO", G3SkyMap::COSMO)
	    .value("none", G3SkyMap::ConvNone) // "None" is reserved in python
	;
	enum_none_converter::from_python<G3SkyMap::MapPolConv, G3SkyMap::ConvNone>();

	bp::class_<G3SkyMap, boost::noncopyable,
	  G3SkyMapPtr>("G3SkyMap",
	  "Base class for 1- and 2-D skymaps of various projections. Usually "
	  "you want a subclass of this (e.g. FlatSkyMap) rather than using it "
	  "directly.", bp::no_init)
	    .def_readonly("__g3frameobject__", true)
	    .def_readwrite("coord_ref", &G3SkyMap::coord_ref,
	      "Coordinate system (maps.MapCoordReference) of the map (e.g. "
	      "Galactic, Equatorial, etc.)")
	    .def_readwrite("pol_type", &G3SkyMap::pol_type,
	      "Polarization type (maps.MapPolType) of the map "
	      "(e.g. maps.MapPolType.Q).")
	    .add_property("pol_conv", &G3SkyMap::GetPolConv, &G3SkyMap::SetPolConv,
	      "Polarization convention (maps.MapPolConv) of the map "
	      "(e.g. maps.MapPolConv.IAU or maps.MapPolConv.COSMO). "
	      "Switching between IAU and COSMO conventions for a U map "
	      "multiplies the U map by -1.")
	    .def_readwrite("units", &G3SkyMap::units,
	      "Unit class (core.G3TimestreamUnits) of the map (e.g. "
	      "core.G3TimestreamUnits.Tcmb). Within each unit class, further "
	      "conversions, for example from K to uK, should use core.G3Units.")
	    .def_readwrite("weighted", &G3SkyMap::weighted,
	      "True if map is multiplied by weights")
	    .add_property("size", &G3SkyMap::size, "Number of pixels in map")
	    .def("__len__", &G3SkyMap::size, "Number of pixels in map")
	    .add_property("shape", &skymap_shape, "Shape of map")
	    .add_property("npix_allocated", &G3SkyMap::NpixAllocated,
	      "Number of pixels in map currently stored in memory")
	    .add_property("npix_nonzero", &G3SkyMap::NpixNonZero,
	      "Number of nonzero pixels in map currently stored in memory")
	    .def_readwrite("overflow", &G3SkyMap::overflow,
	      "Combined value of data processed by "
	      "the map maker but outside of the map area")
	    .def("__getitem__", &skymap_getitem)
	    .def("__setitem__", &skymap_setitem)
	    .def("__copy__", &skymap_copy)
	    .def("__deepcopy__", &skymap_copy)
	    .def("copy", &skymap_copy, "Return a copy of the map object")
	    .def("Clone", &G3SkyMap::Clone,
	      (bp::arg("copy_data")=true),
	       "Return a map of the same type, populated with a copy of the data "
	       "if the argument is true (default), empty otherwise.")
	    .def("IsCompatible", &G3SkyMap::IsCompatible,
	      "Returns true if the input argument is a map with matching dimensions "
	      "and boundaries on the sky.")

	    .def("angles_to_pixels", &G3SkyMap::AnglesToPixels,
	      (bp::arg("alphas"), bp::arg("deltas")),
	       "Compute the 1D pixel location for each of the sky coordinates")
	    .def("pixels_to_angles", &skymap_pixels_to_angles,
	      (bp::arg("pixels")),
	       "Compute the sky coordinates of each of the given 1D pixels")
	    .def("pixel_to_angle", 
	      (std::vector<double> (G3SkyMap::*)(size_t) const) 
	      &G3SkyMap::PixelToAngle, bp::arg("pixel"),
	      "Compute the sky coordinates of the given 1D pixel")
	    .def("angle_to_pixel", &G3SkyMap::AngleToPixel,
	      (bp::arg("alpha"), bp::arg("delta")),
	       "Compute the 1D pixel location of the given sky coordinates.")

	    .def("get_interp_values", &G3SkyMap::GetInterpValues,
	      (bp::arg("alphas"), bp::arg("deltas")),
	       "Return the values at each of the input coordinate locations. "
	       "Computes each value using bilinear interpolation over the "
	       "map pixels.")

	    .def("rebin", &G3SkyMap::Rebin, (bp::arg("scale"), bp::arg("norm")=true),
	      "Rebin the map into larger pixels by summing (if norm is false) "
	      "or averaging (if norm is true) scale-x-scale blocks of pixels "
	      "together.  Returns a new map object.  Map dimensions must be a "
	      "multiple of the rebinning scale.")

	    .def("compact", &G3SkyMap::Compact, (bp::arg("zero_nans")=false),
	      "Convert the map to its default sparse representation, removing "
	      "empty pixels, and optionally also removing NaN values. A map "
	      "that is already sparse will be compactified in place in its "
	      "current representation without additional memory overhead.")

	    .def(bp::self += bp::self)
	    .def(bp::self *= bp::self)
	    .def(bp::self -= bp::self)
	    .def(bp::self /= bp::self)
	    .def(bp::self += double())
	    .def(bp::self *= double())
	    .def(bp::self -= double())
	    .def(bp::self /= double())

	    .def("__add__", &pyskymap_add)
	    .def("__add__", &pyskymap_addd)
	    .def("__radd__", &pyskymap_addd)
	    .def("__sub__", &pyskymap_sub)
	    .def("__sub__", &pyskymap_subd)
	    .def("__rsub__", &pyskymap_rsubd)
	    .def("__mul__", &pyskymap_mult)
	    .def("__mul__", &pyskymap_multd)
	    .def("__rmul__", &pyskymap_multd)
	    .def("__div__", &pyskymap_div)
	    .def("__div__", &pyskymap_divd)
	    .def("__rdiv__", &pyskymap_rdivd)
	    .def("__truediv__", &pyskymap_div)
	    .def("__truediv__", &pyskymap_divd)
	    .def("__rtruediv__", &pyskymap_rdivd)
	    .def("__neg__", &pyskymap_neg)
	;

	EXPORT_FRAMEOBJECT(G3SkyMapWeights, init<>(),
	    "Polarized (Mueller matrix) or unpolarized (scalar) map pixel weights.")
	    .def(bp::init<G3SkyMapConstPtr, bool>(
	      (bp::arg("skymap"), bp::arg("polarized") = true)))
	    .def_readwrite("TT",&G3SkyMapWeights::TT, "Mueller matrix component map")
	    .def_readwrite("TQ",&G3SkyMapWeights::TQ, "Mueller matrix component map")
	    .def_readwrite("TU",&G3SkyMapWeights::TU, "Mueller matrix component map")
	    .def_readwrite("QQ",&G3SkyMapWeights::QQ, "Mueller matrix component map")
	    .def_readwrite("QU",&G3SkyMapWeights::QU, "Mueller matrix component map")
	    .def_readwrite("UU",&G3SkyMapWeights::UU, "Mueller matrix component map")
	    .add_property("size", &G3SkyMapWeights::size, "Number of pixels in weights")
	    .def("__len__", &G3SkyMapWeights::size, "Number of pixels in weights")
	    .def("__copy__", &skymapweights_copy)
	    .def("__deepcopy__", &skymapweights_copy)
	    .def("copy", &skymapweights_copy, "Return a copy of the weights object")
	    .add_property("shape", &skymapweights_shape, "Shape of weights")
	    .add_property("polarized", &G3SkyMapWeights::IsPolarized,
	      "True if all components are set, False if only the TT component is set")
	    .add_property("congruent", &G3SkyMapWeights::IsCongruent,
	      "True if all components are internally compatible with each other")
	    .def("rebin", &G3SkyMapWeights::Rebin, (bp::arg("scale")),
	      "Rebin the weights into larger pixels by summing scale-x-scale blocks "
	      "of pixels together.  Returns a new weights object.  Map dimensions "
	      "must be a multiple of the  rebinning scale.")
	    .def("compact", &G3SkyMapWeights::Compact, (bp::arg("zero_nans")=false),
	      "Convert the map to its default sparse representation, excluding "
	      "empty pixels, and optionally converting NaN values to zeroes.")
	    .def("det", &G3SkyMapWeights::Det,
	      "Return the determinant of the Mueller matrix for each pixel")
	    .def("cond", &G3SkyMapWeights::Cond,
	      "Return the condition number of the Mueller matrix for each pixel")

	    .def("Clone", &G3SkyMapWeights::Clone, (bp::arg("copy_data")=true),
	       "Return weights of the same type, populated with a copy of the data "
	       "if the argument is true (default), empty otherwise.")

	    .def(bp::self += bp::self)
	    .def(bp::self *= FlatSkyMap())
	    .def(bp::self *= HealpixSkyMap())
	    .def(bp::self *= double())
	    .def(bp::self /= double())

	    .def("__add__", &pyskymapweights_add)
	    .def("__mul__", &pyskymapweights_multm)
	    .def("__mul__", &pyskymapweights_multd)
	    .def("__rmul__", &pyskymapweights_multd)
	    .def("__div__", &pyskymapweights_divd)
	    .def("__truediv__", &pyskymapweights_divd)
	;
	register_pointer_conversions<G3SkyMapWeights>();

}

