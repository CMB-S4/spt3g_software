#include <pybindings.h>
#include <serialization.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <G3Logging.h>
#include <G3Units.h>

#include <maps/pointing.h>
#include <maps/FlatSkyProjection.h>

#define COS cos
#define SIN sin
#define ASIN asin
#define ATAN2 atan2

using namespace G3Units;

FlatSkyProjection::FlatSkyProjection(size_t xpix, size_t ypix, double res,
    double alpha_center, double delta_center, double x_res, MapProjection proj,
    double x_center, double y_center)
{
	initialize(xpix, ypix, res, alpha_center, delta_center, x_res, proj,
	    x_center, y_center);
}

FlatSkyProjection::FlatSkyProjection()
{
	initialize(0, 0, 0);
}

FlatSkyProjection::FlatSkyProjection(const FlatSkyProjection & fp)
{
	initialize(fp.xpix_, fp.ypix_, fp.y_res_, fp.alpha0_, fp.delta0_,
	    fp.x_res_, fp.proj_, fp.x0_, fp.y0_);
}

template <class A> void FlatSkyProjection::load(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("xdim", xpix_);
	ar & make_nvp("ydim", ypix_);
	ar & make_nvp("proj", proj_);
	ar & make_nvp("alpha_center", alpha0_);
	ar & make_nvp("delta_center", delta0_);
	// Fix a stupid pickling mistake
	if (v == 1) {
		ar & make_nvp("y_res", y_res_);
		ar & make_nvp("x_res", x_res_);
	} else {
		ar & make_nvp("x_res", x_res_);
		ar & make_nvp("y_res", y_res_);
	}

	if (v > 2) {
		ar & make_nvp("x_center", x0_);
		ar & make_nvp("y_center", y0_);
		if (v == 3) {
			x0_ -= 1;
			y0_ -= 1;
		}
	} else {
		x0_ = 0.0 / 0.0;
		y0_ = 0.0 / 0.0;
	}

	initialize(xpix_, ypix_, y_res_, alpha0_, delta0_, x_res_, proj_, x0_, y0_);
}

template <class A> void FlatSkyProjection::save(A &ar, unsigned v) const
{
	using namespace cereal;

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("xdim", xpix_);
	ar & make_nvp("ydim", ypix_);
	ar & make_nvp("proj", proj_);
	ar & make_nvp("alpha_center", alpha0_);
	ar & make_nvp("delta_center", delta0_);
	ar & make_nvp("x_res", x_res_);
	ar & make_nvp("y_res", y_res_);
	ar & make_nvp("x_center", x0_);
	ar & make_nvp("y_center", y0_);
}

std::string FlatSkyProjection::Description() const
{
	std::ostringstream os;

	os.precision(4);

	os << xpix_ << " x " << ypix_ <<
	    " (" << xpix_ * x_res_ / deg << " x "
	     << ypix_ * y_res_ / deg << " deg) ";

	switch (proj_) {
	case ProjSFL:
		os << "SFL";
		break;
	case ProjCAR:
		os << "CAR";
		break;
	case ProjSIN:
		os << "SIN";
		break;
	case ProjSTG:
		os << "STG";
		break;
	case ProjZEA:
		os << "ZEA";
		break;
	case ProjTAN:
		os << "TAN";
		break;
	case ProjCEA:
		os << "CEA";
		break;
	case ProjBICEP:
		os << "BICEP";
		break;
	default:
		os << "other (" << proj_ << ")";
	}

	os << " centered at (" << x0_ <<
	    ", " << y0_ << ")";

	os << " = (" << alpha0_ / deg << ", "
	    << delta0_ / deg << " deg)";

	return os.str();
}

bool FlatSkyProjection::IsCompatible(const FlatSkyProjection & other) const
{
	bool check = ((xpix_ == other.xpix_) &&
		      (ypix_ == other.ypix_) &&
		      (fabs(x_res_ - other.x_res_) < 1e-8) &&
		      (fabs(y_res_ - other.y_res_) < 1e-8));
	if (proj_ != other.proj_ && (proj_ == ProjNone || other.proj_ == ProjNone)) {
		log_warn("Checking compatibility of maps with projections %d and %d. "
			 "In the future, comparison to a map with projection %d "
			 "(ProjNone) will raise an error.", proj_, other.proj_, ProjNone);
		return check;
	}
	double da = fmod(fabs(alpha0_ - other.alpha0_), 360 * deg);
	if (da > 180 * deg)
		da = 360 * deg - da;
	return (check &&
		(proj_ == other.proj_) &&
		(fabs(delta0_ - other.delta0_) < 1e-8) &&
		(da < 1e-8) &&
		(fabs(x0_ - other.x0_) < 1e-8) &&
		(fabs(y0_ - other.y0_) < 1e-8));
}

void FlatSkyProjection::initialize(size_t xpix, size_t ypix, double res,
    double alpha_center, double delta_center, double x_res, MapProjection proj,
    double x_center, double y_center)
{
	xpix_ = xpix;
	ypix_ = ypix;
	SetProj(proj);
	SetRes(res, x_res);
	SetAngleCenter(alpha_center, delta_center);
	SetXYCenter(x_center, y_center);
}

void FlatSkyProjection::SetProj(MapProjection proj)
{
	proj_ = proj;

	switch (proj_) {
	case Proj0:
	case Proj1:
	case Proj7:
	case Proj9:
		cyl_ = true;
		break;
	default:
		cyl_ = false;
		break;
	}
}

void FlatSkyProjection::SetAlphaCenter(double alpha)
{
	if (alpha < 0)
		alpha += 360 * deg;
	alpha0_ = alpha;
	q0_ = get_origin_rotator(alpha0_, delta0_);
}

void FlatSkyProjection::SetDeltaCenter(double delta)
{
	if (fabs(delta) > 90 * deg)
		log_fatal("Delta center out of range");
	delta0_ = delta;
	sindelta0_ = SIN(delta0_ / rad);
	cosdelta0_ = COS(delta0_ / rad);
	q0_ = get_origin_rotator(alpha0_, delta0_);
}

void FlatSkyProjection::SetAngleCenter(double alpha, double delta)
{
	SetAlphaCenter(alpha);
	SetDeltaCenter(delta);
}

void FlatSkyProjection::SetXCenter(double x)
{
	x0_ = (x != x) ? (xpix_ / 2.0 - 0.5) : x;
}

void FlatSkyProjection::SetYCenter(double y)
{
	y0_ = (y != y) ? (ypix_ / 2.0 - 0.5) : y;
}

void FlatSkyProjection::SetXYCenter(double x, double y)
{
	SetXCenter(x);
	SetYCenter(y);
}

void FlatSkyProjection::SetXRes(double res)
{
	x_res_ = (res != res || res == 0.0) ? y_res_ : res;
}

void FlatSkyProjection::SetYRes(double res)
{
	y_res_ = res;
}

void FlatSkyProjection::SetRes(double res, double x_res)
{
	SetYRes(res);
	SetXRes(x_res);
}

size_t
FlatSkyProjection::XYToPixel(double x, double y) const
{
	// Truncate X/Y coordinates to integer pixels and wrap to 1D.
	// If the pixel is outside of the projection area, return -1.
	// The floor() is important to properly truncate negative numbers,
	// which otherwise pile up at zero from both directions.

	ssize_t ix = (ssize_t)floor(x + 0.5);
	ssize_t iy = (ssize_t)floor(y + 0.5);
	return (ix < 0 || iy < 0 || ix >= (ssize_t) xpix_ || iy >= (ssize_t) ypix_) ?
	    -1 : ix + iy * xpix_;
}

std::vector<double>
FlatSkyProjection::PixelToXY(size_t pixel) const
{
	std::vector<double> out(2, -1);

	if (pixel >= 0 && pixel < xpix_ * ypix_) {
		out[0] = (ssize_t)(pixel % xpix_);
		out[1] = (ssize_t)(pixel / xpix_);
	}

	return out;
}

std::vector<double>
FlatSkyProjection::XYToAngle(double x, double y) const
{
	if (!cyl_) {
		G3Quat q = XYToQuat(x, y);
		double alpha, delta;
		quat_to_ang(q, alpha, delta);
		return {alpha, delta};
	}

	x = (x0_ - x) * x_res_;
	y = (y0_ - y) * y_res_;

	double alpha, delta;

	switch(proj_) {
	case Proj0: {
		delta = delta0_ - y;
		alpha = x / COS(delta / rad) + alpha0_;
		break;
	}
	case Proj1: {
		delta = delta0_ - y;
		alpha = x + alpha0_;
		break;
	}
	case Proj7: {
		delta = ASIN(sindelta0_ - y) * rad;
		alpha = x + alpha0_;
		break;
	}
	case Proj9: {
		delta = delta0_ - y;
		alpha = x / cosdelta0_ + alpha0_;
		break;
	}
	default:
		log_fatal("Proj %d not implemented", proj_);
		break;
	}

	if (alpha < 0)
		alpha += 360 * deg;

	return {alpha, delta};
}

std::vector<double>
FlatSkyProjection::AngleToXY(double alpha, double delta) const
{
	if (!cyl_) {
		G3Quat q = ang_to_quat(alpha, delta);
		return QuatToXY(q);
	}

	static const double circ = 360 * deg;
	static const double halfcirc = 180 * deg;
	static const double quartercirc = 90 * deg;
	double dalpha = alpha - alpha0_;

	if (fabs(delta) > quartercirc)
		return {-1, -1};

	if (dalpha > halfcirc)
		alpha -= circ;
	if (dalpha < -halfcirc)
		alpha += circ;
	dalpha = alpha - alpha0_;

	double x, y;

	switch(proj_) {
	case Proj0: {
		x = dalpha * COS(delta / rad);
		y = delta0_ - delta;
		break;
	}
	case Proj1: {
		x = dalpha;
		y = delta0_ - delta;
		break;
	}
	case Proj7: {
		x = dalpha;
		y = sindelta0_ - SIN(delta / rad);
		break;
	}
	case Proj9: {
		x = dalpha * cosdelta0_;
		y = delta0_ - delta;
		break;
	}
	default:
		log_fatal("Proj %d not implemented", proj_);
		break;
	}

	x = x0_ - x / x_res_;
	y = y0_ - y / y_res_;

	return {x, y};
}

std::vector<double>
FlatSkyProjection::QuatToXY(const G3Quat &q) const
{
	if (cyl_) {
		double a, d;
		quat_to_ang(q, a, d);
		return AngleToXY(a, d);
	}

	// Rotate to projection center
	G3Quat qr = ~q0_ * q * q0_;
	double cc = qr.b();
	double k;

	switch(proj_) {
	case Proj2: {
		k = rad;
		break;
	}
	case Proj3: {
		k = rad * acos(cc) / sqrt((1. - cc) * (1. + cc));
		break;
	}
	case Proj4: {
		k = rad * (2. / (1. + cc));
		break;
	}
	case Proj5: {
		k = rad * sqrt(2. / (1. + cc));
		break;
	}
	case Proj6: {
		k = rad / cc;
		break;
	}
	default:
		log_fatal("Proj %d not implemented", proj_);
		break;
	}

	double x = k * qr.c();
	double y = -k * qr.d();

	x = x0_ - x / x_res_;
	y = y0_ - y / y_res_;

	return {x, y};
}

G3Quat
FlatSkyProjection::XYToQuat(double x, double y) const
{
	if (cyl_) {
		std::vector<double> alphadelta = XYToAngle(x, y);
		return ang_to_quat(alphadelta[0], alphadelta[1]);
	}

	x = (x0_ - x) * x_res_ / rad;
	y = (y0_ - y) * y_res_ / rad;

	double rho = sqrt(x * x + y * y);
	G3Quat q;

	if (rho < 1e-8) {
		q = G3Quat(0, 1, 0, 0);
	} else if (proj_ == Proj2) {
		double cc = sqrt((1. - rho) * (1. + rho));
		q = G3Quat(0, cc, x, -y);
	} else {
		double c;

		switch (proj_) {
		case Proj3:
			c = rho;
			break;
		case Proj4:
			c = 2. * atan(rho / 2.);
			break;
		case Proj5:
			c = 2. * ASIN(rho / 2.);
			break;
		case Proj6:
			c = atan(rho);
			break;
		default:
			log_fatal("Proj %d not implemented", proj_);
			break;
		}

		double cc = COS(c);
		double sc = SIN(c) / rho;
		q = G3Quat(0, cc, x * sc, -y * sc);
	}

	// Rotate from projection center
	return q0_ * q * ~q0_;
}

G3Quat
FlatSkyProjection::PixelToQuat(size_t pixel) const
{
	if (pixel < 0 || pixel >= xpix_ * ypix_)
		return G3Quat(0, 1, 0, 0);
	std::vector<double> xy = PixelToXY(pixel);
	return XYToQuat(xy[0], xy[1]);
}

size_t
FlatSkyProjection::QuatToPixel(const G3Quat &q) const
{
	std::vector<double> xy = QuatToXY(q);
	return XYToPixel(xy[0], xy[1]);
}

std::vector<double>
FlatSkyProjection::PixelToAngle(size_t pixel) const
{
	if (pixel < 0 || pixel >= xpix_ * ypix_)
		return {0., 0.};
	std::vector<double> xy = PixelToXY(pixel);
	return XYToAngle(xy[0], xy[1]);
}

size_t
FlatSkyProjection::AngleToPixel(double alpha, double delta) const
{
	std::vector<double> xy = AngleToXY(alpha, delta);
	return XYToPixel(xy[0], xy[1]);
}

G3VectorQuat
FlatSkyProjection::GetRebinQuats(size_t pixel, size_t scale) const
{
	G3VectorQuat quats(scale * scale, G3Quat(0, 1, 0, 0));

	if (pixel < 0 || pixel >= xpix_ * ypix_) {
		log_debug("Point lies outside of pixel grid\n");
		quats.resize(0);
		return quats;
	}

	std::vector <double> xy = PixelToXY(pixel);
	double x0 = xy[0] - 0.5;
	double y0 = xy[1] - 0.5;

	size_t count = 0;
	for (size_t j = 0; j < scale; j++) {
		double y = y0 + (j + 0.5) / (double) scale;
		for (size_t i = 0; i < scale; i++) {
			double x = x0 + (i + 0.5) / (double) scale;
			quats[count] = XYToQuat(x, y);
			count++;
		}
	}

	return quats;
}

std::vector<double>
FlatSkyProjection::XYToAngleGrad(double x, double y, double h) const
{
	// step size
	double hh = 2 * h;
	const double circ = 360 * deg;
	const double halfcirc = 180 * deg;

	// gradient along x
	std::vector<double> ang1 = XYToAngle(x - h, y);
	std::vector<double> ang2 = XYToAngle(x + h, y);
	if (fabs(ang2[0] - ang1[0]) > halfcirc) {
		ang1[0] = fmod(ang1[0] + halfcirc, circ);
		ang2[0] = fmod(ang2[0] + halfcirc, circ);
	}
	double dax = (ang2[0] - ang1[0]) / hh;
	double ddx = (ang2[1] - ang1[1]) / hh;

	// gradient along y
	std::vector<double> ang3 = XYToAngle(x, y - h);
	std::vector<double> ang4 = XYToAngle(x, y + h);
	if (fabs(ang4[0] - ang3[0]) > halfcirc) {
		ang3[0] = fmod(ang3[0] + halfcirc, circ);
		ang4[0] = fmod(ang4[0] + halfcirc, circ);
	}
	double day = (ang4[0] - ang3[0]) / hh;
	double ddy = (ang4[1] - ang3[1]) / hh;

	return {dax, day, ddx, ddy};
}

std::vector<double>
FlatSkyProjection::PixelToAngleGrad(size_t pixel, double h) const
{
	if (pixel < 0 || pixel >= xpix_ * ypix_)
		return {0., 0., 0., 0.};
	std::vector<double> xy = PixelToXY(pixel);
	return XYToAngleGrad(xy[0], xy[1], h);
}

void FlatSkyProjection::GetInterpPixelsWeights(const G3Quat &q,
    std::vector<size_t> & pixels, std::vector<double> & weights) const
{
	std::vector<double> xy = QuatToXY(q);
	double x = xy[0];
	double y = xy[1];

	pixels = std::vector<size_t>(4, (size_t) -1);
	weights = std::vector<double>(4, 0);

	ssize_t x_1 = (ssize_t)floorf(x);
	ssize_t x_2 = x_1 + 1;
	ssize_t y_1 = (ssize_t)floorf(y);
	ssize_t y_2 = y_1 + 1;
	if (x_1 < 0 || y_1 < 0 || x_2 >= (ssize_t) xpix_ || y_2 >= (ssize_t) ypix_){
		log_debug("Point lies outside of pixel grid\n");
		return;
	}

	pixels[0] = x_1 + y_1 * xpix_;  weights[0] = (x_2 - x) * (y_2 - y);
	pixels[1] = x_2 + y_1 * xpix_;  weights[1] = (x - x_1) * (y_2 - y);
	pixels[2] = x_1 + y_2 * xpix_;  weights[2] = (x_2 - x) * (y - y_1);
	pixels[3] = x_2 + y_2 * xpix_;  weights[3] = (x - x_1) * (y - y_1);
}

std::vector<size_t>
FlatSkyProjection::QueryDisc(const G3Quat &q, double radius) const
{
	static const size_t npts = 72;
	double dd = -2.0 * radius / sqrt(2.0);
	G3Quat qd = get_origin_rotator(0, dd);
	G3Quat p = qd * q * ~qd;
	double pva = q.b();
	double pvb = q.c();
	double pvc = q.d();

	ssize_t xmin = xpix_;
	ssize_t xmax = 0;
	ssize_t ymin = ypix_;
	ssize_t ymax = 0;

	double step = M_PI / (double)(npts);

	// select square patch around center
	for (size_t i = 0; i < npts; i++) {
		double phi = i * step;
		double c = COS(phi);
		double s = SIN(phi);
		G3Quat qv = G3Quat(c, pva * s, pvb * s, pvc * s);
		auto xy = QuatToXY(qv * p * ~qv);
		ssize_t fx = std::floor(xy[0]);
		ssize_t cx = std::ceil(xy[0]);
		ssize_t fy = std::floor(xy[1]);
		ssize_t cy = std::ceil(xy[1]);
		if (fx < xmin)
			xmin = fx < 0 ? 0 : fx;
		if (cx > xmax)
			xmax = cx > (ssize_t) xpix_ ? xpix_ : cx;
		if (fy < ymin)
			ymin = fy < 0 ? 0 : fy;
		if (cy > ymax)
			ymax = cy > (ssize_t) ypix_ ? ypix_ : cy;
	}

	double crad = cos(radius / rad);

	std::vector<size_t> pixels;
	for (ssize_t x = xmin; x < xmax; x++) {
		for (ssize_t y = ymin; y < ymax; y++) {
			size_t pixel = y * xpix_ + x;
			if (pixel > xpix_ * ypix_)
				continue;
			G3Quat qp = PixelToQuat(pixel);
			if (dot3(qp, q) > crad)
				pixels.push_back(pixel);
		}
	}
	std::sort(pixels.begin(), pixels.end());

	return pixels;
}

FlatSkyProjection FlatSkyProjection::Rebin(size_t scale, double x_center, double y_center) const
{
	FlatSkyProjection fp(*this);

	if (scale <= 1)
		return fp;

	fp.xpix_ /= scale;
	fp.ypix_ /= scale;
	fp.SetRes(y_res_ * scale, x_res_ * scale);

	// Use correct center for extracted patch
	if (x_center != x_center && x0_ != (xpix_ / 2.0 - 0.5))
		x_center = (x0_ - xpix_ / 2) / scale + fp.xpix_ / 2;

	if (y_center != y_center && y0_ != (ypix_ / 2.0 - 0.5))
		y_center = (y0_ - ypix_ / 2) / scale + fp.ypix_ / 2;

	fp.SetXYCenter(x_center, y_center);

	return fp;
}

size_t FlatSkyProjection::RebinPixel(size_t pixel, size_t scale) const
{
	size_t x = (pixel % xpix_) / scale;
	size_t y = ((size_t)(pixel / xpix_)) / scale;
	return x + y * (xpix_ / scale);
}

FlatSkyProjection FlatSkyProjection::OverlayPatch(double x0, double y0,
    size_t width, size_t height) const
{
	FlatSkyProjection fp(*this);

	fp.xpix_ = width;
	fp.ypix_ = height;
	fp.SetXYCenter(x0_ - x0 + width / 2, y0_ - y0 + height / 2);

	return fp;
}

std::vector<double> FlatSkyProjection::GetPatchCenter(const FlatSkyProjection &proj) const
{
	// check that input projection is compatible aside from patch location
	FlatSkyProjection fp(proj);
	fp.xpix_ = xpix_;
	fp.ypix_ = ypix_;
	fp.SetXYCenter(x0_, y0_);
	g3_assert(IsCompatible(fp));

	double x0 = x0_ - proj.x0_ + proj.xpix_ / 2;
	double y0 = y0_ - proj.y0_ + proj.ypix_ / 2;

	return {x0, y0};
}

G3_SPLIT_SERIALIZABLE_CODE(FlatSkyProjection);

PYBINDINGS("maps")
{
	using namespace boost::python;

	bp::enum_<MapProjection>("MapProjection")
	    .value("Proj0", Proj0)
	    .value("Proj1", Proj1)
	    .value("Proj2", Proj2)
	    .value("Proj3", Proj3)
	    .value("Proj4", Proj4)
	    .value("Proj5", Proj5)
	    .value("Proj6", Proj6)
	    .value("Proj7", Proj7)
	    .value("Proj8", Proj8)
	    .value("Proj9", Proj9)

	    .value("ProjSansonFlamsteed", ProjSansonFlamsteed)
	    .value("ProjSFL", ProjSFL)
	    .value("ProjPlateCarree", ProjPlateCarree)
	    .value("ProjCAR", ProjCAR)
	    .value("ProjOrthographic", ProjOrthographic)
	    .value("ProjSIN", ProjSIN)
	    .value("ProjStereographic", ProjStereographic)
	    .value("ProjSTG", ProjSTG)
	    .value("ProjLambertAzimuthalEqualArea",
	      ProjLambertAzimuthalEqualArea)
	    .value("ProjZEA", ProjZEA)
	    .value("ProjGnomonic", ProjGnomonic)
	    .value("ProjTAN", ProjTAN)
	    .value("ProjCylindricalEqualArea",
	      ProjCylindricalEqualArea)
	    .value("ProjCEA", ProjCEA)
	    .value("ProjBICEP", ProjBICEP)
	    .value("ProjNone", ProjNone)
	;
	enum_none_converter::from_python<MapProjection, ProjNone>();
}
