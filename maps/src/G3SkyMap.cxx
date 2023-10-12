#include <pybindings.h>
#include <serialization.h>
#include <G3Units.h>
#include <maps/G3SkyMap.h>
#include <maps/G3SkyMapMask.h>
#include <maps/FlatSkyMap.h>
#include <maps/HealpixSkyMap.h>
#include <maps/pointing.h>
#include <cmath>

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

G3SkyMapPtr
G3SkyMap::ArrayClone(boost::python::object val) const
{

	G3SkyMapPtr skymap = Clone(false);
	skymap->FillFromArray(val);
	return skymap;
}


G3SkyMapWeights::G3SkyMapWeights(G3SkyMapConstPtr ref, bool polarized) :
    TT(ref->Clone(false)), TQ(polarized ? ref->Clone(false) : NULL),
    TU(polarized ? ref->Clone(false) : NULL),
    QQ(polarized ? ref->Clone(false) : NULL),
    QU(polarized ? ref->Clone(false) : NULL),
    UU(polarized ? ref->Clone(false) : NULL)
{
	TT->pol_type = G3SkyMap::None;
	if (polarized){
		TQ->pol_type = G3SkyMap::None;
		TU->pol_type = G3SkyMap::None;
		QQ->pol_type = G3SkyMap::None;
		QU->pol_type = G3SkyMap::None;
		UU->pol_type = G3SkyMap::None;
	}
}

G3SkyMapWeights::G3SkyMapWeights(const G3SkyMapWeights &r, bool copy_data) :
    TT(r.TT->Clone(copy_data)), TQ(!r.TQ ? NULL : r.TQ->Clone(copy_data)),
    TU(!r.TU ? NULL : r.TU->Clone(copy_data)),
    QQ(!r.QQ ? NULL : r.QQ->Clone(copy_data)),
    QU(!r.QU ? NULL : r.QU->Clone(copy_data)),
    UU(!r.UU ? NULL : r.UU->Clone(copy_data))
{
}

std::vector<size_t>
G3SkyMap::AnglesToPixels(const std::vector<double> & alphas,
    const std::vector<double> & deltas) const
{
	std::vector<size_t> pixels(alphas.size());

	for (size_t i = 0; i < alphas.size(); i++) {
		pixels[i] = AngleToPixel(alphas[i], deltas[i]);
	}

	return pixels;
}

void
G3SkyMap::PixelsToAngles(const std::vector<size_t> & pixels,
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

size_t
G3SkyMap::AngleToPixel(double alpha, double delta) const
{
	quat q = ang_to_quat(alpha, delta);

	return QuatToPixel(q);
}

std::vector<double>
G3SkyMap::PixelToAngle(size_t pixel) const
{
	quat q = PixelToQuat(pixel);
	double alpha, delta;
	quat_to_ang(q, alpha, delta);

	return {alpha, delta};
}

std::vector<size_t>
G3SkyMap::QuatsToPixels(const G3VectorQuat &quats) const
{
	std::vector<size_t> pixels(quats.size());
	for (size_t i = 0; i < quats.size(); i++)
		pixels[i] = QuatToPixel(quats[i]);

	return pixels;
}

G3VectorQuat
G3SkyMap::PixelsToQuats(const std::vector<size_t> &pixels) const
{
	G3VectorQuat quats(pixels.size());
	for (size_t i = 0; i < pixels.size(); i++)
		quats[i] = PixelToQuat(pixels[i]);

	return quats;
}


static boost::python::tuple
skymap_pixels_to_angles(const G3SkyMap & skymap,
    const std::vector<size_t> & pixels)
{
	std::vector<double> alphas, deltas;
	skymap.PixelsToAngles(pixels, alphas, deltas);

	return boost::python::make_tuple(alphas, deltas);
}

static boost::python::tuple
skymap_pixel_to_angle(const G3SkyMap & skymap, size_t pixel)
{
	std::vector<double> alphadelta = skymap.PixelToAngle(pixel);

	return boost::python::make_tuple(alphadelta[0], alphadelta[1]);
}


size_t G3SkyMap::size() const
{
	size_t s = 1;
	for (size_t i : shape())
		s *= i;
	return s;
}

G3SkyMapMaskPtr
G3SkyMap::MakeMask(bool zero_nans, bool zero_infs) const
{
	G3SkyMapMaskPtr m(new G3SkyMapMask(*this, true, zero_nans, zero_infs));
	return m;
}

G3SkyMap &G3SkyMap::operator+=(const G3SkyMap & rhs)
{
	g3_assert(IsCompatible(rhs));
	g3_assert(units == rhs.units);
	g3_assert(weighted == rhs.weighted);

	for (size_t i = 0; i < rhs.size(); i++)
		(*this)[i] += rhs.at(i);
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
	g3_assert(units == rhs.units);
	g3_assert(weighted == rhs.weighted);

	for (size_t i = 0; i < rhs.size(); i++)
		(*this)[i] -= rhs.at(i);
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
	if (units == G3Timestream::None)
		units = rhs.units;
	if (rhs.weighted and !(weighted))
		weighted = true;

	for (size_t i = 0; i < rhs.size(); i++)
		(*this)[i] *= rhs.at(i);
	return *this;
}

G3SkyMap &G3SkyMap::operator*=(const G3SkyMapMask &rhs)
{
	g3_assert(rhs.IsCompatible(*this));

	for (size_t i = 0; i < size(); i++) {
		if (!rhs.at(i) && this->at(i) != 0)
			(*this)[i] = 0;
	}
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
	if (units == G3Timestream::None)
		units = rhs.units;
	if (rhs.weighted and !(weighted))
		weighted = true;

	for (size_t i = 0; i < rhs.size(); i++)
		(*this)[i] /= rhs.at(i);
	return *this;
}

G3SkyMap &G3SkyMap::operator/=(double rhs)
{
	for (size_t i = 0; i < size(); i++)
		(*this)[i] /= rhs;
	return *this;
}


void
G3SkyMap::GetInterpPixelsWeights(double alpha, double delta,
    std::vector<size_t> & pixels, std::vector<double> & weights) const
{
	quat q = ang_to_quat(alpha, delta);
	GetInterpPixelsWeights(q, pixels, weights);
}

void
G3SkyMap::GetRebinAngles(size_t pixel, size_t scale,
    std::vector<double> &alphas, std::vector<double> &deltas) const
{
	auto quats = GetRebinQuats(pixel, scale);
	alphas = std::vector<double>(quats.size());
	deltas = std::vector<double>(quats.size());

	for (size_t i = 0; i < quats.size(); i++) {
		double alpha, delta;
		quat_to_ang(quats[i], alpha, delta);
		alphas[i] = alpha;
		deltas[i] = delta;
	}
}

double
G3SkyMap::GetInterpPrecalc(const std::vector<size_t> & pix,
    const std::vector<double> & weight) const
{
	double outval = 0;
	for (size_t i = 0; i < pix.size(); i++) {
		outval += this->at(pix[i]) * weight[i];
	}
	return outval;
}

double
G3SkyMap::GetInterpValue(double alpha, double delta) const
{
	quat q = ang_to_quat(alpha, delta);
	return GetInterpValue(q);
}

std::vector<double>
G3SkyMap::GetInterpValues(const std::vector<double> & alphas,
    const std::vector<double> & deltas) const
{
	std::vector<double> outvals(alphas.size());

	for (size_t i = 0; i < alphas.size(); i++) {
		outvals[i] = GetInterpValue(alphas[i], deltas[i]);
	}

	return outvals;
}

double
G3SkyMap::GetInterpValue(quat q) const
{
	std::vector<size_t> pix;
	std::vector<double> weight;
	GetInterpPixelsWeights(q, pix, weight);
	return GetInterpPrecalc(pix, weight);
}

std::vector<double>
G3SkyMap::GetInterpValues(const G3VectorQuat & quats) const
{
	std::vector<double> outvals(quats.size());

	for (size_t i = 0; i < quats.size(); i++)
		outvals[i] = GetInterpValue(quats[i]);

	return outvals;
}

std::vector<size_t>
G3SkyMap::QueryDisc(double alpha, double delta, double radius) const
{
	quat q = ang_to_quat(alpha, delta);
	return QueryDisc(q, radius);
}

std::vector<size_t>
G3SkyMap::QueryAlphaEllipse(double alpha ,double delta, double a, double b) const
{
	quat q = ang_to_quat(alpha, delta);
	return QueryAlphaEllipse(q, a, b);
}

std::vector<size_t>
G3SkyMap::QueryAlphaEllipse(quat q, double a, double b) const
{
	double rmaj = a > b ? a : b;
	double rmin = a > b ? b : a;
	double sd = q.R_component_4();
	double cd = sqrt((1 - sd) * (1 + sd));
	double da = ACOS(COS(rmaj) / COS(rmin)) / cd;

	quat qda = get_origin_rotator(da, 0);
	quat ql = qda * q * ~qda;
	quat qr = ~qda * q * qda;

	auto disc = QueryDisc(q, rmaj);

	std::vector<size_t> pixels;
	for (auto i: disc) {
		quat qp = PixelToQuat(i);
		double d = ACOS(dot3(qp, ql)) + ACOS(dot3(qp, qr));
		if (d < 2 * rmaj)
			pixels.push_back(i);
	}

	return pixels;
}

static double
skymap_getitem(const G3SkyMap &skymap, ssize_t i)
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
skymap_setitem(G3SkyMap &skymap, ssize_t i, double val)
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

	if (pol_conv != pol_conv_) {
		log_warn("Sign of U map flipped in changing pol_conv from %s to %s",
		    pol_conv == G3SkyMap::IAU ? "COSMO" : "IAU",
		    pol_conv == G3SkyMap::IAU ? "IAU" : "COSMO");
		(*this) *= -1;
        }

	pol_conv_ = pol_conv;
}

static bp::tuple
skymap_shape(const G3SkyMap &skymap)
{
	// Swap to match numpy's convention for shape()
	std::vector<size_t> shape = skymap.shape();
	std::vector<uint64_t> pyshape;
	for (ssize_t i = shape.size() - 1; i >= 0; i--)
		pyshape.push_back(shape[i]);
	
	return bp::tuple(pyshape);
}

static bp::tuple
skymapweights_shape(const G3SkyMapWeights &weights)
{
	return skymap_shape(*(weights.TT));
}

static G3SkyMapPtr
skymap_copy(const G3SkyMap &r)
{
	return r.Clone(true);
}

static G3SkyMapWeightsPtr
skymapweights_copy(const G3SkyMapWeights &r)
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
pyskymap_multm(const G3SkyMap &a, const G3SkyMapMask &b)
{
	g3_assert(b.IsCompatible(a));
	G3SkyMapPtr rv = a.Clone(false);
	for (auto i: b) {
		if (b.at(i.first) && a.at(i.first) != 0)
			(*rv)[i.first] = a.at(i.first);
	}

	return rv;
}

static G3SkyMapPtr
pyskymap_imultm(G3SkyMapPtr a, const G3SkyMapMask &b)
{
	(*a) *= b;
	return a;
}

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

static void
pyskymap_ipowd(G3SkyMap &a, double b)
{
	if (b == 0) {
		a *= 0;
		a += 1;
	} else {
		for (size_t i = 0; i < a.size(); i++) {
			double v = a.at(i);
			if (v != 0)
				a[i] = pow(v, b);
		}
	}
}

static G3SkyMapPtr
pyskymap_powd(const G3SkyMap &a, double b)
{
	G3SkyMapPtr rv;
	if (b == 0) {
		rv = a.Clone(false);
		(*rv) += 1;
	} else {
		rv = a.Clone(true);
		pyskymap_ipowd(*rv, b);
	}
	return rv;
}

static void
pyskymap_ipow(G3SkyMap &a, const G3SkyMap &b)
{
	g3_assert(a.IsCompatible(b));
	g3_assert(b.units == G3Timestream::None);
	for (size_t i = 0; i < a.size(); i++) {
		double va = a.at(i);
		double vb = b.at(i);
		if (va != 0 || vb == 0) {
			a[i] = pow(va, vb);
		}
	}
}

static G3SkyMapPtr
pyskymap_pow(const G3SkyMap &a, const G3SkyMap &b)
{
	G3SkyMapPtr rv = a.Clone(true);
	pyskymap_ipow(*rv, b);
	return rv;
}

#define skymap_comp(oper) \
G3SkyMapMask \
G3SkyMap::operator oper(const G3SkyMap &rhs) \
{ \
	g3_assert(IsCompatible(rhs)); \
	g3_assert(units == rhs.units); \
	G3SkyMapMask rv(*this); \
	for (size_t i = 0; i < size(); i++) { \
		if (at(i) oper rhs.at(i)) \
			rv[i] = true; \
	} \
	return rv; \
}

skymap_comp(<)
skymap_comp(<=)
skymap_comp(==)
skymap_comp(!=)
skymap_comp(>=)
skymap_comp(>)

#define skymap_compd(oper) \
G3SkyMapMask \
G3SkyMap::operator oper(double rhs) \
{ \
	G3SkyMapMask rv(*this);	\
	for (size_t i = 0; i < size(); i++) { \
		if (at(i) oper rhs) \
			rv[i] = true; \
	} \
	return rv; \
}

skymap_compd(<)
skymap_compd(<=)
skymap_compd(==)
skymap_compd(!=)
skymap_compd(>=)
skymap_compd(>)

bool
G3SkyMap::all(G3SkyMapMaskConstPtr where) const
{
	if (!where) {
		for (size_t i = 0; i < size(); i++) {
			if (at(i) == 0)
				return false;
		}
	} else {
		g3_assert(where->IsCompatible(*this));

		for (size_t i = 0; i < size(); i++) {
			if (!where->at(i))
				continue;
			if (at(i) == 0)
				return false;
		}
	}

	return true;
}

bool
G3SkyMap::any(G3SkyMapMaskConstPtr where) const
{
	if (!where) {
		for (size_t i = 0; i < size(); i++) {
			if (at(i) != 0)
				return true;
		}
	} else {
		g3_assert(where->IsCompatible(*this));

		for (size_t i = 0; i < size(); i++) {
			if (!where->at(i))
				continue;
			if (at(i) != 0)
				return true;
		}
	}

	return false;
}

double
G3SkyMap::sum(G3SkyMapMaskConstPtr where) const
{
	double s = 0;

	if (!where) {
		for (size_t i = 0; i < size(); i++) {
			s += at(i);
		}
	} else {
		g3_assert(where->IsCompatible(*this));

		for (size_t i = 0; i < size(); i++) {
			if (!where->at(i))
				continue;
			s += at(i);
		}
	}

	return s;
}

double
G3SkyMap::mean(G3SkyMapMaskConstPtr where) const
{
	double m = 0;
	size_t n = 0;

	if (!where) {
		n = size();
		for (size_t i = 0; i < n; i++) {
			m += at(i);
		}
	} else {
		g3_assert(where->IsCompatible(*this));

		for (size_t i = 0; i < size(); i++) {
			if (!where->at(i))
				continue;
			m += at(i);
			n++;
		}
	}

	return m / n;
}

double
G3SkyMap::var(size_t ddof, G3SkyMapMaskConstPtr where) const
{
	double m1 = 0;
	double m2 = 0;
	double a, b;
	size_t n = 0;

	if (!where) {
		n = size();
		for (size_t i = 0; i < n; i++) {
			a = at(i);
			m1 += a;
			m2 += a * a;
		}
	} else {
		g3_assert(where->IsCompatible(*this));

		for (size_t i = 0; i < size(); i++) {
			if (!where->at(i))
				continue;
			n++;
			a = at(i);
			m1 += a;
			m2 += a * a;
		}
	}

	return (m2 - m1 * m1 / (double)n) / (double)(n - ddof);
}

double
G3SkyMap::std(size_t ddof, G3SkyMapMaskConstPtr where) const
{
	return sqrt(var(ddof, where));
}

double
G3SkyMap::min(G3SkyMapMaskConstPtr where) const
{
	double m = std::numeric_limits<double>::infinity();
	double v;

	if (!where) {
		for (size_t i = 0; i < size(); i++) {
			v = at(i);
			if (v < m)
				m = v;
		}
	} else {
		g3_assert(where->IsCompatible(*this));

		for (size_t i = 0; i < size(); i++) {
			if (!where->at(i))
				continue;
			v = at(i);
			if (v < m)
				m = v;
		}
	}

	return m;
}

size_t
G3SkyMap::argmin(G3SkyMapMaskConstPtr where) const
{
	double m = std::numeric_limits<double>::infinity();
	double v;
	size_t j = 0;

	if (!where) {
		for (size_t i = 0; i < size(); i++) {
			v = at(i);
			if (v < m) {
				m = v;
				j = i;
			}
		}
	} else {
		g3_assert(where->IsCompatible(*this));

		for (size_t i = 0; i < size(); i++) {
			if (!where->at(i))
				continue;
			v = at(i);
			if (v < m) {
				m = v;
				j = i;
			}
		}
	}

	return j;
}

double
G3SkyMap::max(G3SkyMapMaskConstPtr where) const
{
	double m = -1 * std::numeric_limits<double>::infinity();
	double v;

	if (!where) {
		for (size_t i = 0; i < size(); i++) {
			v = at(i);
			if (v > m)
				m = v;
		}
	} else {
		g3_assert(where->IsCompatible(*this));

		for (size_t i = 0; i < size(); i++) {
			if (!where->at(i))
				continue;
			v = at(i);
			if (v > m)
				m = v;
		}
	}

	return m;
}

size_t
G3SkyMap::argmax(G3SkyMapMaskConstPtr where) const
{
	double m = -1 * std::numeric_limits<double>::infinity();
	double v;
	size_t j = 0;

	if (!where) {
		for (size_t i = 0; i < size(); i++) {
			v = at(i);
			if (v > m) {
				m = v;
				j = i;
			}
		}
	} else {
		g3_assert(where->IsCompatible(*this));

		for (size_t i = 0; i < size(); i++) {
			if (!where->at(i))
				continue;
			v = at(i);
			if (v > m) {
				m = v;
				j = i;
			}
		}
	}

	return j;
}

double
G3SkyMap::median(G3SkyMapMaskConstPtr where) const
{
	std::vector<double> data;
	size_t npix = where ? where->sum() : size();
	if (npix == 0)
		return 0;

	if (!where) {
		for (size_t i = 0; i < size(); i++) {
			data.push_back(at(i));
		}
	} else {
		g3_assert(where->IsCompatible(*this));

		for (size_t i = 0; i < size(); i++) {
			if (!where->at(i))
				continue;
			data.push_back(at(i));
		}
	}

	npix = data.size();
	size_t n = npix / 2;

	std::nth_element(data.begin(), data.begin() + n, data.end());

	// odd-length array
	if (npix % 2)
		return data[n];

	// even-length array
	double m = data[n];
	std::nth_element(data.begin(), data.begin() + n - 1, data.end());
	return (m + data[n - 1]) / 2.;
}

#define skymap_test(oper) \
G3SkyMapMask \
G3SkyMap::oper(G3SkyMapMaskConstPtr where) const \
{ \
	G3SkyMapMask mask(*this); \
	if (!where) { \
		for (size_t i = 0; i < size(); i++) { \
			if (std::oper(at(i))) \
				mask[i] = true; \
		} \
	} else { \
		g3_assert(where->IsCompatible(*this)); \
		for (size_t i = 0; i < size(); i++) { \
			if (!where->at(i)) \
				continue; \
			if (std::oper(at(i))) \
				mask[i] = true; \
		} \
	} \
	return mask; \
}

skymap_test(isinf);
skymap_test(isnan);
skymap_test(isfinite);

#define skymap_nanoper(oper, type) \
type \
G3SkyMap::nan##oper(G3SkyMapMaskConstPtr where) const \
{ \
	G3SkyMapMask mask = isnan(where); \
	mask.invert(); \
	return oper(boost::make_shared<G3SkyMapMask>(mask)); \
}

skymap_nanoper(sum, double);
skymap_nanoper(mean, double);
skymap_nanoper(median, double);
skymap_nanoper(min, double);
skymap_nanoper(max, double);
skymap_nanoper(argmin, size_t);
skymap_nanoper(argmax, size_t);

double
G3SkyMap::nanvar(size_t ddof, G3SkyMapMaskConstPtr where) const
{
	G3SkyMapMask mask = isnan(where);
	mask.invert();
	return var(ddof, boost::make_shared<G3SkyMapMask>(mask));
}

double
G3SkyMap::nanstd(size_t ddof, G3SkyMapMaskConstPtr where) const
{
	return sqrt(nanvar(ddof, where));
}

static bool
pyskymap_bool(G3SkyMap &skymap)
{
	PyErr_SetString(PyExc_ValueError,
	    "ValueError: The truth value of a G3SkyMap is ambiguous. Use m.any() or m.all()");
	bp::throw_error_already_set();

	return false;
}

static std::vector<uint64_t>
pyskymap_nonzero(G3SkyMap &skymap)
{
	return skymap.MakeMask()->NonZeroPixels();
}

void
G3SkyMap::ApplyMask(const G3SkyMapMask &mask, bool inverse)
{
	g3_assert(mask.IsCompatible(*this));

	for (size_t i = 0; i < size(); i++) {
		if (!at(i))
			continue;
		if (mask.at(i) == inverse)
			(*this)[i] = 0;
	}
}

void
G3SkyMapWeights::ApplyMask(const G3SkyMapMask &mask, bool inverse)
{
	if (!!TT)
		TT->ApplyMask(mask, inverse);
	if (!!TQ)
		TQ->ApplyMask(mask, inverse);
	if (!!TU)
		TU->ApplyMask(mask, inverse);
	if (!!QQ)
		QQ->ApplyMask(mask, inverse);
	if (!!QU)
		QU->ApplyMask(mask, inverse);
	if (!!UU)
		UU->ApplyMask(mask, inverse);
}

G3SkyMapWeights &
G3SkyMapWeights::operator+=(const G3SkyMapWeights &rhs)
{
	g3_assert(IsPolarized() == rhs.IsPolarized());

	if (!!TT)
		*TT += *rhs.TT;
	if (!!TQ)
		*TQ += *rhs.TQ;
	if (!!TU)
		*TU += *rhs.TU;
	if (!!QQ)
		*QQ += *rhs.QQ;
	if (!!QU)
		*QU += *rhs.QU;
	if (!!UU)
		*UU += *rhs.UU;

	return *this;
}

#define skymapweights_inplace(op, rhs_type) \
G3SkyMapWeights &G3SkyMapWeights::operator op(rhs_type rhs) \
{ \
	if (!!TT) \
		*TT op rhs; \
	if (!!TQ) \
		*TQ op rhs; \
	if (!!TU) \
		*TU op rhs; \
	if (!!QQ) \
		*QQ op rhs; \
	if (!!QU) \
		*QU op rhs; \
	if (!!UU) \
		*UU op rhs; \
	return *this; \
}

skymapweights_inplace(*=, const G3SkyMap &);
skymapweights_inplace(*=, const G3SkyMapMask &);
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

static G3SkyMapWeightsPtr
pyskymapweights_multma(const G3SkyMapWeights &a, const G3SkyMapMask &b)
{
	G3SkyMapWeightsPtr rv(new G3SkyMapWeights());
	if (!!a.TT)
		rv->TT = pyskymap_multm(*a.TT, b);
	if (!!a.TQ)
		rv->TQ = pyskymap_multm(*a.TQ, b);
	if (!!a.TU)
		rv->TU = pyskymap_multm(*a.TU, b);
	if (!!a.QQ)
		rv->QQ = pyskymap_multm(*a.QQ, b);
	if (!!a.QU)
		rv->QU = pyskymap_multm(*a.QU, b);
	if (!!a.UU)
		rv->UU = pyskymap_multm(*a.UU, b);

	return rv;
}

static G3SkyMapWeightsPtr
pyskymapweights_imultma(G3SkyMapWeightsPtr a, const G3SkyMapMask &b)
{
	if (!!a->TT)
		(*a->TT) *= b;
	if (!!a->TQ)
		(*a->TQ) *= b;
	if (!!a->TU)
		(*a->TU) *= b;
	if (!!a->QQ)
		(*a->QQ) *= b;
	if (!!a->QU)
		(*a->QU) *= b;
	if (!!a->UU)
		(*a->UU) *= b;

	return a;
}

MuellerMatrix MuellerMatrix::Inv() const
{
	MuellerMatrix m;
	double c = Cond();
	if (tt == 0 || c != c || c > 1e12) {
		if ((c != c || c > 1e12) && tt != 0)
			log_trace("Singular matrix found when inverting!  Cond is %lE\n", c);
		m.tt = m.tq = m.tu = m.qq = m.qu = m.uu = 0.0 / 0.0;
		return m;
	}

	double d = Det();
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
	if (p1 == 0) {
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

	MuellerMatrix B;
	B.tt = ttq;
	B.qq = qqq;
	B.uu = uuq;
	B.tq = tq;
	B.tu = tu;
	B.qu = qu;
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
	double c = lmax / lmin;
	if (c < 0)
		return 0.0 / 0.0;
	return c;
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

G3SkyMapWeightsPtr G3SkyMapWeights::Inv() const
{
	G3SkyMapWeightsPtr out = this->Clone(false);
	out->TT->ConvertToDense();
	if (!!TQ)
		out->TQ->ConvertToDense();
	if (!!TU)
		out->TU->ConvertToDense();
	if (!!QQ)
		out->QQ->ConvertToDense();
	if (!!QU)
		out->QU->ConvertToDense();
	if (!!UU)
		out->UU->ConvertToDense();

	for (size_t i = 0; i < TT->size(); i++)
		(*out)[i] = this->at(i).Inv();

	return out;
}

G3SkyMapWeightsPtr G3SkyMapWeights::Rebin(size_t scale) const
{
	g3_assert(IsCongruent());
	G3SkyMapWeightsPtr out(new G3SkyMapWeights());

	out->TT = !TT ? NULL : TT->Rebin(scale, false);
	out->TQ = !TQ ? NULL : TQ->Rebin(scale, false);
	out->TU = !TU ? NULL : TU->Rebin(scale, false);
	out->QQ = !QQ ? NULL : QQ->Rebin(scale, false);
	out->QU = !QU ? NULL : QU->Rebin(scale, false);
	out->UU = !UU ? NULL : UU->Rebin(scale, false);

	return out;
}

void G3SkyMapWeights::Compact(bool zero_nans)
{
	g3_assert(IsCongruent());

	if (!!TT)
		TT->Compact(zero_nans);
	if (!!TQ)
		TQ->Compact(zero_nans);
	if (!!TU)
		TU->Compact(zero_nans);
	if (!!QQ)
		QQ->Compact(zero_nans);
	if (!!QU)
		QU->Compact(zero_nans);
	if (!!UU)
		UU->Compact(zero_nans);
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
	    .value("E", G3SkyMap::E)
	    .value("B", G3SkyMap::B)
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
	    .def("copy", &skymap_copy, "Return a copy of the map object")
	    .def("clone", &G3SkyMap::Clone,
	      (bp::arg("copy_data")=true),
	       "Return a map of the same type, populated with a copy of the data "
	       "if the argument is true (default), empty otherwise.")
	    .def("array_clone", &G3SkyMap::ArrayClone,
	      (bp::arg("array")),
	       "Return a map of the same type, populated with a copy of the input "
	       "numpy array")
	    .def("compatible", &G3SkyMap::IsCompatible,
	      "Returns true if the input argument is a map with matching dimensions "
	      "and boundaries on the sky.")
	    .def("nonzero", &pyskymap_nonzero, "Return indices of non-zero pixels in the map")

	    .def("angles_to_pixels", &G3SkyMap::AnglesToPixels,
	      (bp::arg("alphas"), bp::arg("deltas")),
	       "Compute the 1D pixel location for each of the sky coordinates "
	       "(vectorized)")
	    .def("pixels_to_angles", &skymap_pixels_to_angles,
	      (bp::arg("pixels")),
	       "Compute the sky coordinates of each of the given 1D pixels "
	       "(vectorized)")
	    .def("angle_to_pixel", &G3SkyMap::AnglesToPixels,
	      (bp::arg("alphas"), bp::arg("deltas")),
	       "Compute the 1D pixel location for each of the sky coordinates "
	       "(vectorized)")
	    .def("pixel_to_angle", &skymap_pixels_to_angles,
	      (bp::arg("pixels")),
	       "Compute the sky coordinates of each of the given 1D pixels "
	       "(vectorized)")
	    .def("quats_to_pixels", &G3SkyMap::QuatsToPixels,
	       "Compute the 1D pixel location for each of the sky coordinates "
	       "(vectorized), expressed as quaternion rotations from the pole.")
	    .def("pixels_to_quats", &G3SkyMap::PixelsToQuats,
	       "Compute the sky coordinates, expressed as quaternion rotations "
	       "from the pole, for each of the given 1-D pixel coordinates.")
	    .def("angle_to_pixel", &G3SkyMap::AngleToPixel,
	      (bp::arg("alpha"), bp::arg("delta")),
	       "Compute the 1D pixel location of the given sky position.")
	    .def("pixel_to_angle", &skymap_pixel_to_angle, (bp::arg("pixel")),
	       "Compute the sky coordinates (alpha, delta) of the given 1D pixel")
	    .def("quat_to_pixel", &G3SkyMap::QuatToPixel,
	       "Compute the 1D pixel location of the given sky position, "
	       "expressed as a quaternion rotation from the pole.")
	    .def("pixel_to_quat", &G3SkyMap::PixelToQuat,
	       "Compute the quaternion rotation from the pole corresponding to "
	       "the given 1D pixel.")
	    .def("quat_to_pixel", &G3SkyMap::QuatsToPixels,
	       "Compute the 1D pixel location for each of the sky coordinates "
	       "(vectorized), expressed as quaternion rotations from the pole.")
	    .def("pixel_to_quat", &G3SkyMap::PixelsToQuats,
	       "Compute the sky coordinates, expressed as quaternion rotations "
	       "from the pole, for each of the given 1-D pixel coordinates.")

	    .def("query_disc",
	      (std::vector<size_t> (G3SkyMap::*)(double, double, double) const)
		&G3SkyMap::QueryDisc,
	       (bp::arg("alpha"), bp::arg("delta"), bp::arg("radius")),
	       "Return a list of pixel indices whose centers are located within "
	       "a disc of the given radius at the given sky coordinates.")

	    .def("query_disc",
	      (std::vector<size_t> (G3SkyMap::*)(quat, double) const)
		&G3SkyMap::QueryDisc,
	       (bp::arg("quat"), bp::arg("radius")),
	       "Return a list of pixel indices whose centers are located within "
	       "a disc of the given radius at the given sky coordinates.")

	    .def("query_alpha_ellipse",
	      (std::vector<size_t> (G3SkyMap::*)(double, double, double, double) const)
		&G3SkyMap::QueryAlphaEllipse,
	       (bp::arg("alpha"), bp::arg("delta"), bp::arg("a"), bp::arg("b")),
	       "Return a list of pixel indices whose centers are located within an "
	       "ellipse extended in the alpha direction, at the given alpha and "
	       "delta sky coordinates, with semimajor and semiminor axes a and b.")

	    .def("query_alpha_ellipse",
	      (std::vector<size_t> (G3SkyMap::*)(quat, double, double) const)
		&G3SkyMap::QueryAlphaEllipse,
	       (bp::arg("quat"), bp::arg("a"), bp::arg("b")),
	       "Return a list of pixel indices whose centers are located within an "
	       "ellipse extended in the alpha direction, at the given alpha and "
	       "delta sky coordinates, with semimajor and semiminor axes a and b.")

	    .def("get_interp_values",
	      (std::vector<double> (G3SkyMap::*)(const std::vector<double> &,
		const std::vector<double> &) const) &G3SkyMap::GetInterpValues,
	      (bp::arg("alphas"), bp::arg("deltas")),
	       "Return the values at each of the input coordinate locations. "
	       "Computes each value using bilinear interpolation over the "
	       "map pixels.")

	    .def("get_interp_values",
	      (std::vector<double> (G3SkyMap::*)(const G3VectorQuat &) const)
	      &G3SkyMap::GetInterpValues, (bp::arg("quats")),
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

	    .def("to_mask", &G3SkyMap::MakeMask,
	      (bp::arg("zero_nans")=false, bp::arg("zero_infs")=false),
	      "Create a G3SkyMapMask object from the parent map, with pixels "
	      "set to true where the map is non-zero (and optionally non-nan "
	      "and/or finite).")

	    .def("apply_mask", &G3SkyMap::ApplyMask,
	      (bp::arg("mask"), bp::arg("inverse")=false),
	      "Apply a mask in-place to the map, optionally inverting which "
	      "pixels are zeroed.  If inverse = False, this is equivalent to "
	      "in-place multiplication by the mask.")

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
	    .def("__imul__", &pyskymap_imultm)
	    .def("__mul__", &pyskymap_mult)
	    .def("__mul__", &pyskymap_multm)
	    .def("__mul__", &pyskymap_multd)
	    .def("__rmul__", &pyskymap_multm)
	    .def("__rmul__", &pyskymap_multd)
	    .def("__div__", &pyskymap_div)
	    .def("__div__", &pyskymap_divd)
	    .def("__rdiv__", &pyskymap_rdivd)
	    .def("__truediv__", &pyskymap_div)
	    .def("__truediv__", &pyskymap_divd)
	    .def("__rtruediv__", &pyskymap_rdivd)
	    .def("__neg__", &pyskymap_neg)
	    .def("__pos__", &skymap_copy)
	    .def("__pow__", &pyskymap_pow)
	    .def("__pow__", &pyskymap_powd)

	    .def(bp::self < bp::self)
	    .def(bp::self <= bp::self)
	    .def(bp::self == bp::self)
	    .def(bp::self != bp::self)
	    .def(bp::self >= bp::self)
	    .def(bp::self > bp::self)
	    .def(bp::self < double())
	    .def(bp::self <= double())
	    .def(bp::self == double())
	    .def(bp::self != double())
	    .def(bp::self >= double())
	    .def(bp::self > double())

	    .def("__bool__", &pyskymap_bool)
	    .def("_cany", &G3SkyMap::any, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("_call", &G3SkyMap::all, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("_csum", &G3SkyMap::sum, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("_cmean", &G3SkyMap::mean, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("median", &G3SkyMap::median, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("_cvar", &G3SkyMap::var, (bp::arg("ddof")=0,
	        bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("_cstd", &G3SkyMap::std, (bp::arg("ddof")=0,
	        bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("_cmin", &G3SkyMap::min, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("_cmax", &G3SkyMap::max, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("_cargmin", &G3SkyMap::argmin, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("_cargmax", &G3SkyMap::argmax, (bp::arg("where")=G3SkyMapMaskConstPtr()))

	    .def("nansum", &G3SkyMap::nansum, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("nanmean", &G3SkyMap::nanmean, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("nanmedian", &G3SkyMap::nanmedian, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("nanvar", &G3SkyMap::nanvar, (bp::arg("ddof")=0,
	        bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("nanstd", &G3SkyMap::nanstd, (bp::arg("ddof")=0,
	        bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("nanmin", &G3SkyMap::nanmin, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("nanmax", &G3SkyMap::nanmax, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("nanargmin", &G3SkyMap::nanargmin, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("nanargmax", &G3SkyMap::nanargmax, (bp::arg("where")=G3SkyMapMaskConstPtr()))

	    .def("isinf", &G3SkyMap::isinf, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("isnan", &G3SkyMap::isnan, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	    .def("isfinite", &G3SkyMap::isfinite, (bp::arg("where")=G3SkyMapMaskConstPtr()))
	;
	boost::python::implicitly_convertible<G3SkyMapPtr, G3SkyMapConstPtr>();

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
	    .def("copy", &skymapweights_copy, "Return a copy of the weights object")
	    .add_property("shape", &skymapweights_shape, "Shape of weights")
	    .add_property("polarized", &G3SkyMapWeights::IsPolarized,
	      "True if all components are set, False if only the TT component is set")
	    .add_property("congruent", &G3SkyMapWeights::IsCongruent,
	      "True if all components are internally compatible with each other")
	    .def("compatible", &G3SkyMapWeights::IsCompatible,
	      "Returns true if the input argument is a map with matching dimensions "
	      "and boundaries on the sky.")
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
	    .def("inv", &G3SkyMapWeights::Inv,
	      "Return the inverse of the Mueller matrix for each pixel")

	    .def("apply_mask", &G3SkyMapWeights::ApplyMask,
	      (bp::arg("mask"), bp::arg("inverse")=false),
	      "Apply a mask in-place to the weights, optionally inverting which "
	      "pixels are zeroed.  If inverse = False, this is equivalent to "
	      "in-place multiplication by the mask.")

	    .def("clone", &G3SkyMapWeights::Clone, (bp::arg("copy_data")=true),
	       "Return weights of the same type, populated with a copy of the data "
	       "if the argument is true (default), empty otherwise.")

	    .def(bp::self += bp::self)
	    .def(bp::self *= FlatSkyMap())
	    .def(bp::self *= HealpixSkyMap())
	    .def(bp::self *= double())
	    .def(bp::self /= double())

	    .def("__add__", &pyskymapweights_add)
	    .def("__imul__", &pyskymapweights_imultma)
	    .def("__mul__", &pyskymapweights_multm)
	    .def("__mul__", &pyskymapweights_multma)
	    .def("__mul__", &pyskymapweights_multd)
	    .def("__rmul__", &pyskymapweights_multm)
	    .def("__rmul__", &pyskymapweights_multma)
	    .def("__rmul__", &pyskymapweights_multd)
	    .def("__div__", &pyskymapweights_divd)
	    .def("__truediv__", &pyskymapweights_divd)
	;
	register_pointer_conversions<G3SkyMapWeights>();

}

