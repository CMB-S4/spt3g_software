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
G3_SERIALIZABLE_CODE(G3SkyMapWithWeights);

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

G3SkyMapWithWeights::G3SkyMapWithWeights(G3SkyMapConstPtr ref,
  bool weighted, bool polarized, std::string map_id_) :
    T(ref->Clone(false)), Q(polarized ? ref->Clone(false) : NULL),
    U(polarized ? ref->Clone(false) : NULL),
    weights(weighted ? new G3SkyMapWeights(ref, polarized) : NULL),
    map_id(map_id_)
{
}

G3SkyMapWithWeights::G3SkyMapWithWeights(const G3SkyMapWithWeights &r) :
    T(r.T->Clone(true)), Q(!r.Q ? NULL : r.Q->Clone(true)),
    U(!r.U ? NULL : r.U->Clone(true)),
    weights(!r.weights ? NULL : r.weights->Clone(true)),
    map_id(r.map_id)
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

G3SkyMapWithWeights &
G3SkyMapWithWeights::operator+=(const G3SkyMapWithWeights &rhs)
{
	g3_assert(IsPolarized() == rhs.IsPolarized());
	g3_assert(IsWeighted() == rhs.IsWeighted());

	*T += *rhs.T;

	if (IsPolarized()) {
		*Q += *rhs.Q;
		*U += *rhs.U;
	}

	if (IsWeighted())
		*weights += *rhs.weights;

	return *this;
}

#define skymapwithweights_inplace(op, rhs_type) \
G3SkyMapWithWeights &G3SkyMapWithWeights::operator op(rhs_type rhs) \
{ \
	*T op rhs; \
	if (IsPolarized()) { \
		*Q op rhs; \
		*U op rhs; \
	} \
	if (IsWeighted()) \
		*weights op rhs; \
	return *this; \
}

skymapwithweights_inplace(*=, const G3SkyMap &);
skymapwithweights_inplace(*=, double);
skymapwithweights_inplace(/=, double);

#define skymapwithweights_pynoninplace(name, oper, rhs_type) \
static G3SkyMapWithWeightsPtr \
pyskymapwithweights_##name(const G3SkyMapWithWeights &a, const rhs_type b) \
{ \
	G3SkyMapWithWeightsPtr rv = a.Clone(true); \
	(*rv) oper b; \
	return rv; \
}

skymapwithweights_pynoninplace(add, +=, G3SkyMapWithWeights &);
skymapwithweights_pynoninplace(multm, *=, G3SkyMap &);
skymapwithweights_pynoninplace(multd, *=, double);
skymapwithweights_pynoninplace(divd, /=, double);

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
	double p = (tt - q)*(tt - q) + (qq - q)*(qq - q) + (uu - q)*(uu - q) + 2.*p1;
	p = sqrt(p / 6.0);

	MuellerMatrix Q;
	Q.tt = Q.qq = Q.uu = q;
	MuellerMatrix B = *this;
	B -= Q;
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

void G3SkyMapWithWeights::ApplyWeights(G3SkyMapWeightsPtr w)
{
	g3_assert(!IsWeighted());
	g3_assert(IsCongruent());
	g3_assert(w->IsCongruent());
	g3_assert(T->IsCompatible(*(w->TT)));

	if (IsPolarized()) {
		for (size_t pix = 0; pix < T->size(); pix++) {
			StokesVector v = this->at(pix);
			if (!(v.t == 0 && v.q == 0 && v.u == 0))
				(*this)[pix] = w->at(pix) * v;
		}
	} else {
		(*T) *= *(w->TT);
	}

	T->weighted = true;
	if (IsPolarized()) {
		Q->weighted = true;
		U->weighted = true;
	}

	// Store pointer to weights here
	weights = w;
}

G3SkyMapWeightsPtr G3SkyMapWithWeights::RemoveWeights()
{
	g3_assert(IsWeighted());
	g3_assert(IsCongruent());

	// Dividing by empty weights fills with nans
	T->ConvertToDense();
	if (IsPolarized()) {
		Q->ConvertToDense();
		U->ConvertToDense();
		for (size_t pix = 0; pix < T->size(); pix++)
			(*this)[pix] /= weights->at(pix);
	} else {
		(*T) /= *(weights->TT);
	}

	T->weighted = false;
	if (IsPolarized()) {
		Q->weighted = false;
		U->weighted = false;
	}

	// Remove pointer to weights
	G3SkyMapWeightsPtr wout = weights;
	weights.reset();

	return wout;
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

G3SkyMapWithWeightsPtr G3SkyMapWithWeights::Rebin(size_t scale) const
{
	g3_assert(IsCongruent());
	G3SkyMapWithWeightsPtr out(new G3SkyMapWithWeights());

	bool w = IsWeighted();
	bool p = IsPolarized();
	out->T = T->Rebin(scale, !w);
	out->Q = p ? Q->Rebin(scale, !w) : NULL;
	out->U = p ? U->Rebin(scale, !w) : NULL;
	out->weights = w ? weights->Rebin(scale) : NULL;
	out->map_id = map_id;

	return out;
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

	bp::class_<G3SkyMap, boost::noncopyable,
	  G3SkyMapPtr>("G3SkyMap",
	  "Base class for 1- and 2-D skymaps of various projections. Usually "
	  "you want a subclass of this (e.g. FlatSkyMap) rather than using it "
	  "directly.", bp::no_init)
	    .def_readwrite("coord_ref", &G3SkyMap::coord_ref,
	      "Coordinate system (maps.MapCoordReference) of the map (e.g. "
	      "Galactic, Equatorial, etc.)")
	    .def_readwrite("pol_type", &G3SkyMap::pol_type,
	      "Polarization type (maps.MapPolType) of the map "
	      "(e.g. maps.MapPolType.Q)")
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

	EXPORT_FRAMEOBJECT(G3SkyMapWeights, init<>(), "generic sky weight")
	    .def(bp::init<G3SkyMapConstPtr, bool>(
	      (bp::arg("skymap"), bp::arg("polarized") = true)))
	    .def_readwrite("TT",&G3SkyMapWeights::TT)
	    .def_readwrite("TQ",&G3SkyMapWeights::TQ)
	    .def_readwrite("TU",&G3SkyMapWeights::TU)
	    .def_readwrite("QQ",&G3SkyMapWeights::QQ)
	    .def_readwrite("QU",&G3SkyMapWeights::QU)
	    .def_readwrite("UU",&G3SkyMapWeights::UU)
	    .add_property("polarized", &G3SkyMapWeights::IsPolarized)
	    .add_property("congruent", &G3SkyMapWeights::IsCongruent)
	    .def("rebin", &G3SkyMapWeights::Rebin)
	    .def("det", &G3SkyMapWeights::Det)
	    .def("cond", &G3SkyMapWeights::Cond)
	    
	    .def("Clone", &G3SkyMapWeights::Clone)

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

	EXPORT_FRAMEOBJECT(G3SkyMapWithWeights, init<>(), "Container for (potentially) polarized maps and weights")
	    .def(bp::init<G3SkyMapPtr, bool, bool, std::string>(
	      (bp::arg("stub_map"),
	       bp::arg("weighted"),
	       bp::arg("polarized"),
	       bp::arg("map_id") = "")))
	    .def_readwrite("T",&G3SkyMapWithWeights::T)
	    .def_readwrite("Q",&G3SkyMapWithWeights::Q)
	    .def_readwrite("U",&G3SkyMapWithWeights::U)
	    .def_readwrite("weights",&G3SkyMapWithWeights::weights)
	    .def_readwrite("map_id",&G3SkyMapWithWeights::map_id)
	    .add_property("weighted",&G3SkyMapWithWeights::IsWeighted)
	    .add_property("polarized",&G3SkyMapWithWeights::IsPolarized)
	    .add_property("congruent",&G3SkyMapWithWeights::IsCongruent)
	    .def("remove_weights",&G3SkyMapWithWeights::RemoveWeights)
	    .def("apply_weights",&G3SkyMapWithWeights::ApplyWeights)
	    .def("rebin", &G3SkyMapWithWeights::Rebin)
	    
	    .def("Clone", &G3SkyMapWithWeights::Clone)

	    .def(bp::self += bp::self)
	    .def(bp::self *= FlatSkyMap())
	    .def(bp::self *= HealpixSkyMap())
	    .def(bp::self *= double())
	    .def(bp::self /= double())

	    .def("__add__", &pyskymapwithweights_add)
	    .def("__mul__", &pyskymapwithweights_multm)
	    .def("__mul__", &pyskymapwithweights_multd)
	    .def("__rmul__", &pyskymapwithweights_multd)
	    .def("__div__", &pyskymapwithweights_divd)
	    .def("__truediv__", &pyskymapwithweights_divd)
	;
	register_pointer_conversions<G3SkyMapWithWeights>();

}

