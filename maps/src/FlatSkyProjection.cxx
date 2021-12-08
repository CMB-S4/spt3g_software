#include <pybindings.h>
#include <serialization.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <G3Logging.h>
#include <G3Units.h>

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
	initialize(fp.xpix_, fp.ypix_, fp.y_res_, fp.alphac_, fp.deltac_,
	    fp.x_res_, fp.proj_, fp.xc_, fp.yc_);
}

template <class A> void FlatSkyProjection::load(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("xdim", xpix_);
	ar & make_nvp("ydim", ypix_);
	ar & make_nvp("proj", proj_);
	ar & make_nvp("alpha_center", alphac_);
	ar & make_nvp("delta_center", deltac_);
	// Fix a stupid pickling mistake
	if (v == 1) {
		ar & make_nvp("y_res", y_res_);
		ar & make_nvp("x_res", x_res_);
	} else {
		ar & make_nvp("x_res", x_res_);
		ar & make_nvp("y_res", y_res_);
	}

	if (v > 2) {
		ar & make_nvp("x_center", xc_);
		ar & make_nvp("y_center", yc_);
		if (v == 3) {
			xc_ -= 1;
			yc_ -= 1;
		}
	} else {
		xc_ = 0.0 / 0.0;
		yc_ = 0.0 / 0.0;
	}

	initialize(xpix_, ypix_, y_res_, alphac_, deltac_, x_res_, proj_, xc_, yc_);
}

template <class A> void FlatSkyProjection::save(A &ar, unsigned v) const
{
	using namespace cereal;

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("xdim", xpix_);
	ar & make_nvp("ydim", ypix_);
	ar & make_nvp("proj", proj_);
	ar & make_nvp("alpha_center", alphac_);
	ar & make_nvp("delta_center", deltac_);
	ar & make_nvp("x_res", x_res_);
	ar & make_nvp("y_res", y_res_);
	ar & make_nvp("x_center", xc_);
	ar & make_nvp("y_center", yc_);
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

	os << " centered at (" << xc_ <<
	    ", " << yc_ << ")";

	os << " = (" << alphac_ / deg << ", "
	    << deltac_ / deg << " deg)";

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
	return (check &&
		(proj_ == other.proj_) &&
		(fabs(delta0_ - other.delta0_) < 1e-8) &&
		(fmod(fabs(alpha0_ - other.alpha0_), 360 * deg) < 1e-8) &&
		(fabs(x0_ - other.x0_) < 1e-8) &&
		(fabs(y0_ - other.y0_) < 1e-8));
}

void FlatSkyProjection::initialize(size_t xpix, size_t ypix, double res,
    double alpha_center, double delta_center, double x_res, MapProjection proj,
    double x_center, double y_center)
{
	init_ = false;
	xpix_ = xpix;
	ypix_ = ypix;
	proj_ = proj;
	SetRes(res, x_res);
	SetAngleCenter(alpha_center, delta_center);
	SetXYCenter(x_center, y_center);
	init_ = true;
}

void FlatSkyProjection::SetProj(MapProjection proj)
{
	proj_ = proj;
	SetAngleCenter(alphac_, deltac_);
	SetXYCenter(xc_, yc_);
}

void FlatSkyProjection::SetAlphaCenter(double alpha)
{
	alphac_ = alpha;
	alpha0_ = alpha;
}

void FlatSkyProjection::SetDeltaCenter(double delta)
{
	if (fabs(delta) > 90 * deg)
		log_fatal("Delta center out of range");
	deltac_ = delta;
	delta0_ = delta;
	if (proj_ == Proj0 || proj_ == Proj1 || proj_ == Proj7 || proj_ == Proj9)
		delta0_ = 0.0;
	sindelta0_ = SIN(delta0_ / rad);
	cosdelta0_ = COS(delta0_ / rad);
	if (proj_ == Proj9) {
		if (fabs(90 * deg - fabs(deltac_)) < 1e-12)
			log_fatal("Projection %d is not valid at the poles", proj_);
		pc_.resize(1);
		pc_[0] = 1. / COS(deltac_ / rad);
	} else {
		pc_.clear();
	}
	if (init_)
		SetYCenter(yc_);
}

void FlatSkyProjection::SetAngleCenter(double alpha, double delta)
{
	SetAlphaCenter(alpha);
	SetDeltaCenter(delta);
}

void FlatSkyProjection::SetXCenter(double x)
{
	xc_ = (x != x) ? (xpix_ / 2.0 - 0.5) : x;
	x0_ = xc_;
}

void FlatSkyProjection::SetYCenter(double y)
{
	yc_ = (y != y) ? (ypix_ / 2.0 - 0.5) : y;
	y0_ = yc_;
	if (proj_ == Proj0 || proj_ == Proj1 || proj_ == Proj7 || proj_ == Proj9)
		y0_ -= (proj_ == Proj7 ? SIN(deltac_ / rad) : deltac_) / y_res_;
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

long
FlatSkyProjection::XYToPixel(double x, double y) const
{
	// Truncate X/Y coordinates to integer pixels and wrap to 1D.
	// If the pixel is outside of the projection area, return -1.
	// The floor() is important to properly truncate negative numbers,
	// which otherwise pile up at zero from both directions.

	long ix = (long)floor(x + 0.5);
	long iy = (long)floor(y + 0.5);
	return (ix < 0 || iy < 0 || ix >= xpix_ || iy >= ypix_) ?
	    -1 : ix + iy * xpix_;
}

std::vector<double>
FlatSkyProjection::PixelToXY(long pixel) const
{
	std::vector<double> out(2, -1);

	if (pixel >= 0 && pixel < xpix_ * ypix_) {
		out[0] = (int)(pixel % xpix_);
		out[1] = (int)(pixel / xpix_);
	}

	return out;
}

std::vector<double>
FlatSkyProjection::XYToAngle(double x, double y, bool wrap_alpha) const
{
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
	case Proj2: {
		x /= rad;
		y /= rad;
		double rho = sqrt(x * x + y * y);
		double a, b;
		if (rho < 1e-8) {
			a = 0;
			b = delta0_;
		} else {
			double c = ASIN(rho);
			double cc = COS(c);
			a = ATAN2(x, cc * cosdelta0_ + y * sindelta0_) * rad;
			b = ASIN(cc * sindelta0_ - y * cosdelta0_) * rad;
		}
		delta = b;
		alpha = a + alpha0_;
		break;
	}
	case Proj3: {
		x /= rad;
		y /= rad;
		double rho = sqrt(x * x + y * y);
		double a, b;
		if (rho < 1e-8) {
			a = 0;
			b = delta0_;
		} else {
			double c = rho;
			double cc = COS(c);
			double sc = SIN(c);
			a = ATAN2(x * sc, rho * cc * cosdelta0_ + y * sc * sindelta0_) * rad;
			b = ASIN(cc * sindelta0_ - y * sc / rho * cosdelta0_) * rad;
		}
		delta = b;
		alpha = a + alpha0_;
		break;
	}
	case Proj4: {
		x /= rad;
		y /= rad;
		double rho = sqrt(x * x + y * y);
		double a, b;
		if (rho < 1e-8) {
			a = 0;
			b = delta0_;
		} else {
			double c = 2. * atan(rho / 2.);
			double cc = COS(c);
			double sc = SIN(c);
			a = ATAN2(x * sc, rho * cc * cosdelta0_ + y * sc * sindelta0_) * rad;
			b = ASIN(cc * sindelta0_ - y * sc / rho * cosdelta0_) * rad;
		}
		delta = b;
		alpha = a + alpha0_;
		break;
	}
	case Proj5: {
		x /= rad;
		y /= rad;
		double rho = sqrt(x * x + y * y);
		double a, b;
		if (rho < 1e-8) {
			a = 0;
			b = delta0_;
		} else {
			double c = 2. * ASIN(rho / 2.);
			double cc = COS(c);
			double sc = SIN(c);
			a = ATAN2(x * sc, rho * cc * cosdelta0_ + y * sc * sindelta0_) * rad;
			b = ASIN(cc * sindelta0_ - y * sc / rho * cosdelta0_) * rad;
		}
		delta = b;
		alpha = a + alpha0_;
		break;
	}
	case Proj6: {
		x /= rad;
		y /= rad;
		double rho = sqrt(x * x + y * y);
		double a, b;
		if (rho < 1e-8) {
			a = 0;
			b = delta0_;
		} else {
			double c = atan(rho);
			double cc = COS(c);
			double sc = SIN(c);
			a = ATAN2(x * sc, rho * cc * cosdelta0_ + y * sc * sindelta0_) * rad;
			b = ASIN(cc * sindelta0_ - y * sc / rho * cosdelta0_) * rad;
		}
		delta = b;
		alpha = a + alpha0_;
		break;
	}
	case Proj7: {
		delta = delta0_ - ASIN(y) * rad;
		alpha = x + alpha0_;
		break;
	}
	case Proj9: {
		delta = delta0_ - y;
		alpha = x * pc_[0] + alpha0_;
		break;
	}
	default:
		log_fatal("Proj %d not implemented", proj_);
		break;
	}

	static const double circ = 360 * deg;
	static const double halfcirc = 180 * deg;
	double dalpha = alpha - alpha0_;

	if (wrap_alpha) {
		alpha = fmod((alpha < 0) ? (circ * (1. + ceilf(fabs(alpha) / circ)) + alpha) : alpha, circ);
	} else {
		if (dalpha > halfcirc)
			alpha -= circ;
		if (dalpha < -halfcirc)
			alpha += circ;
	}

	return {alpha, delta};
}

std::vector<double>
FlatSkyProjection::AngleToXY(double alpha, double delta) const
{
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
	case Proj2: {
		dalpha /= rad;
		delta /= rad;
		double cosdelta = COS(delta);
		y = (cosdelta * sindelta0_ * COS(dalpha) - cosdelta0_ * SIN(delta)) * rad;
		x = cosdelta * SIN(dalpha) * rad;
		break;
	}
	case Proj3: {
		dalpha /= rad;
		delta /= rad;
		double cosdelta = COS(delta);
		double sindelta = SIN(delta);
		double cosdalpha = COS(dalpha);
		double c = ASIN(sindelta0_ * sindelta + cosdelta0_ * cosdelta * cosdalpha);
		double k = (90 * deg - rad * c) / COS(c);
		y = k * (cosdelta * sindelta0_ * cosdalpha - cosdelta0_ * sindelta);
		x = k * cosdelta * SIN(dalpha);
		break;
	}
	case Proj4: {
		dalpha /= rad;
		delta /= rad;
		double cosdelta = COS(delta);
		double sindelta = SIN(delta);
		double cosdalpha = COS(dalpha);
		double k = rad * (2. / (1. + sindelta0_ * sindelta + cosdelta0_ * cosdelta * cosdalpha));
		y = k * (cosdelta * sindelta0_ * cosdalpha - cosdelta0_ * sindelta);
		x = k * cosdelta * SIN(dalpha);
		break;
	}
	case Proj5: {
		dalpha /= rad;
		delta /= rad;
		double cosdelta = COS(delta);
		double sindelta = SIN(delta);
		double cosdalpha = COS(dalpha);
		double k = rad * sqrt(2. / (1. + sindelta0_ * sindelta + cosdelta0_ * cosdelta * cosdalpha));
		y = k * (cosdelta * sindelta0_ * cosdalpha - cosdelta0_ * sindelta);
		x = k * cosdelta * SIN(dalpha);
		break;
	}
	case Proj6: {
		dalpha /= rad;
		delta /= rad;
		double cosdelta = COS(delta);
		double sindelta = SIN(delta);
		double cosdalpha = COS(dalpha);
		double k = rad / (sindelta0_ * sindelta + cosdelta0_ * cosdelta * cosdalpha);
		y = k * (cosdelta * sindelta0_ * cosdalpha - cosdelta0_ * sindelta);
		x = k * cosdelta * SIN(dalpha);
		break;
	}
	case Proj7: {
		x = dalpha;
		y = SIN((delta0_ - delta) / rad);
		break;
	}
	case Proj9: {
		x = dalpha / pc_[0];
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
FlatSkyProjection::PixelToAngle(long pixel, bool wrap_alpha) const
{
	if (pixel < 0 || pixel >= xpix_ * ypix_)
		return {0., 0.};
	std::vector<double> xy = PixelToXY(pixel);
	return XYToAngle(xy[0], xy[1], wrap_alpha);
}

long
FlatSkyProjection::AngleToPixel(double alpha, double delta) const
{
	std::vector<double> xy = AngleToXY(alpha, delta);
	return XYToPixel(xy[0], xy[1]);
}

void FlatSkyProjection::GetRebinAngles(long pixel, size_t scale,
    std::vector<double> & alphas, std::vector<double> & deltas,
    bool wrap_alpha) const
{
	alphas = std::vector<double>(scale * scale);
	deltas = std::vector<double>(scale * scale);

	if (pixel < 0 || pixel >= xpix_ * ypix_) {
		log_debug("Point lies outside of pixel grid\n");
		return;
	}

	std::vector <double> xy = PixelToXY(pixel);
	double x0 = xy[0] - 0.5;
	double y0 = xy[1] - 0.5;

	size_t count = 0;
	for (size_t j = 0; j < scale; j++) {
		double y = y0 + (j + 0.5) / (double) scale;
		for (size_t i = 0; i < scale; i++) {
			double x = x0 + (i + 0.5) / (double) scale;
			std::vector<double> ang = XYToAngle(x, y, wrap_alpha);
			alphas[count] = ang[0];
			deltas[count] = ang[1];
			count++;
		}
	}
}

std::vector<double>
FlatSkyProjection::XYToAngleGrad(double x, double y, double h) const
{
	// step size
	double hh = 2 * h;
	const double circ = 360 * deg;
	const double halfcirc = 180 * deg;

	// gradient along x
	std::vector<double> ang1 = XYToAngle(x - h, y, false);
	std::vector<double> ang2 = XYToAngle(x + h, y, false);
	if (fabs(ang2[0] - ang1[0]) > halfcirc) {
		ang1[0] = fmod(ang1[0] + halfcirc, circ);
		ang2[0] = fmod(ang2[0] + halfcirc, circ);
	}
	double dax = (ang2[0] - ang1[0]) / hh;
	double ddx = (ang2[1] - ang1[1]) / hh;

	// gradient along y
	std::vector<double> ang3 = XYToAngle(x, y - h, false);
	std::vector<double> ang4 = XYToAngle(x, y + h, false);
	if (fabs(ang4[0] - ang3[0]) > halfcirc) {
		ang3[0] = fmod(ang3[0] + halfcirc, circ);
		ang4[0] = fmod(ang4[0] + halfcirc, circ);
	}
	double day = (ang4[0] - ang3[0]) / hh;
	double ddy = (ang4[1] - ang3[1]) / hh;

	return {dax, day, ddx, ddy};
}

std::vector<double>
FlatSkyProjection::PixelToAngleGrad(long pixel, double h) const
{
	if (pixel < 0 || pixel >= xpix_ * ypix_)
		return {0., 0., 0., 0.};
	std::vector<double> xy = PixelToXY(pixel);
	return XYToAngleGrad(xy[0], xy[1], h);
}

void FlatSkyProjection::GetInterpPixelsWeights(double alpha, double delta,
    std::vector<long> & pixels, std::vector<double> & weights) const
{
	std::vector<double> xy = AngleToXY(alpha, delta);
	double x = xy[0];
	double y = xy[1];

	pixels = std::vector<long>(4, -1);
	weights = std::vector<double>(4, 0);

	int x_1 = (int)floorf(x);
	int x_2 = x_1 + 1;
	int y_1 = (int)floorf(y);
	int y_2 = y_1 + 1;
	if (x_1 < 0 || y_1 < 0 || x_2 >= xpix_ || y_2 >= ypix_){
		log_debug("Point lies outside of pixel grid\n");
		return;
	}

	pixels[0] = x_1 + y_1 * xpix_;  weights[0] = (x_2 - x) * (y_2 - y);
	pixels[1] = x_2 + y_1 * xpix_;  weights[1] = (x - x_1) * (y_2 - y);
	pixels[2] = x_1 + y_2 * xpix_;  weights[2] = (x_2 - x) * (y - y_1);
	pixels[3] = x_2 + y_2 * xpix_;  weights[3] = (x - x_1) * (y - y_1);
}

FlatSkyProjection FlatSkyProjection::Rebin(size_t scale, double x_center, double y_center) const
{
	FlatSkyProjection fp(*this);

	if (scale <= 1)
		return fp;

	fp.xpix_ /= scale;
	fp.ypix_ /= scale;
	fp.SetRes(y_res_ * scale, x_res_ * scale);
	fp.SetXYCenter(x_center, y_center);

	return fp;
}

FlatSkyProjection FlatSkyProjection::OverlayPatch(double x0, double y0,
    size_t width, size_t height) const
{
	FlatSkyProjection fp(*this);

	fp.xpix_ = width;
	fp.ypix_ = height;
	fp.SetXYCenter(xc_ - x0 + width / 2, yc_ - y0 + height / 2);

	return fp;
}

std::vector<double> FlatSkyProjection::GetPatchCenter(const FlatSkyProjection &proj) const
{
	// check that input projection is compatible aside from patch location
	FlatSkyProjection fp(proj);
	fp.xpix_ = xpix_;
	fp.ypix_ = ypix_;
	fp.SetXYCenter(xc_, yc_);
	g3_assert(IsCompatible(fp));

	double x0 = xc_ - proj.xc_ + proj.xpix_ / 2;
	double y0 = yc_ - proj.yc_ + proj.ypix_ / 2;

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
