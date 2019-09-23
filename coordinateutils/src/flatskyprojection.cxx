#include <pybindings.h>
#include <serialization.h>

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <G3Logging.h>
#include <G3Units.h>

#include <coordinateutils/flatskyprojection.h>

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

	os.precision(1);

	os << xpix_ << " x " << ypix_ <<
	    " (" << xpix_ * x_res_ / deg << " x "
	     << ypix_ * y_res_ / deg << " deg) ";

	switch (proj_) {
	case ProjSansonFlamsteed:
		os << "Sanson-Flamsteed";
		break;
	case ProjCAR:
		os << "Plate Carree";
		break;
	case ProjSIN:
		os << "Orthographic";
		break;
	case ProjStereographic:
		os << "Stereographic";
		break;
	case ProjLambertAzimuthalEqualArea:
		os << "Lambert Azimuthal Equal Area";
		break;
	case ProjGnomonic:
		os << "Gnomonic";
		break;
	case ProjCEA:
		os << "Cylindrical Equal Area";
		break;
	case ProjBICEP:
		os << "BICEP";
		break;
	default:
		os << "other (" << proj_ << ")";
	}

	return os.str();
}

bool FlatSkyProjection::IsCompatible(const FlatSkyProjection & other) const
{
	return ((xpix_ == other.xpix_) &&
		(ypix_ == other.ypix_) &&
		(proj_ == other.proj_) &&
		(fabs(alpha0_ - other.alpha0_) < 1e-12) &&
		(fabs(delta0_ - other.delta0_) < 1e-12) &&
		(fabs(x_res_ - other.x_res_) < 1e-12) &&
		(fabs(y_res_ - other.y_res_) < 1e-12) &&
		(fabs(x0_ - other.x0_) < 1e-12) &&
		(fabs(y0_ - other.y0_) < 1e-12));
}

void FlatSkyProjection::initialize(size_t xpix, size_t ypix, double res,
    double alpha_center, double delta_center, double x_res, MapProjection proj,
    double x_center, double y_center)
{
	xpix_ = xpix;
	ypix_ = ypix;
	set_proj(proj);
	set_angle_center(alpha_center, delta_center);
	set_res(res, x_res);
	set_xy_center(x_center, y_center);
}

void FlatSkyProjection::set_proj(MapProjection proj)
{
	proj_ = proj;
}

void FlatSkyProjection::set_alpha_center(double alpha)
{
	alpha0_ = alpha;
}

void FlatSkyProjection::set_delta_center(double delta)
{
	delta0_ = delta;
	sindelta0_ = SIN(delta0_ / rad);
	cosdelta0_ = COS(delta0_ / rad);
	if (proj_ == Proj7) {
		delta0_ = fabs(delta0_);
		pv_.resize(1);
		pv_[0] = cosdelta0_ * cosdelta0_ / rad;
	}
}

void FlatSkyProjection::set_angle_center(double alpha, double delta)
{
	set_alpha_center(alpha);
	set_delta_center(delta);
}

void FlatSkyProjection::set_x_center(double x)
{
	x0_ = (x != x) ? (xpix_ / 2.0 + 0.5) : x;
	x_min_ = (x0_ - xpix_) * x_res_;
}

void FlatSkyProjection::set_y_center(double y)
{
	y0_ = (y != y) ? (ypix_ / 2.0 + 0.5) : y;
	if (proj_ == Proj0 || proj_ == Proj1 || proj_ == Proj9) {
		y0_ += delta0_ / y_res_;
		set_delta_center(0.0);
	}
	y_min_ = (y0_ - ypix_) * y_res_;
}

void FlatSkyProjection::set_xy_center(double x, double y)
{
	set_x_center(x);
	set_y_center(y);
}

void FlatSkyProjection::set_xres(double res)
{
	if (res <= 0)
		res = (proj_ == Proj9) ? (y_res_ / cosdelta0_) : y_res_;
	x_res_ = res;
}

void FlatSkyProjection::set_yres(double res)
{
	y_res_ = res;
}

void FlatSkyProjection::set_res(double res, double x_res)
{
	set_yres(res);
	set_xres(x_res);
}

long
FlatSkyProjection::xy_to_pixel(double x, double y) const
{
	long ix = (long) (x + 0.5);
	long iy = (long) (y + 0.5);
	return (ix < 0 || iy < 0 || ix >= xpix_ || iy >= ypix_) ? xpix_ * ypix_ : ix + iy * xpix_;
}

std::vector<double>
FlatSkyProjection::pixel_to_xy(long pixel) const
{
	std::vector<double> out(2, -1);

	if (pixel >= 0 && pixel < xpix_ * ypix_) {
		out[0] = (int)(pixel % xpix_);
		out[1] = (int)(pixel / xpix_);
	}

	return out;
}

std::vector<double>
FlatSkyProjection::xy_to_angle(double x, double y, bool wrap_alpha) const
{
	x = -1. * (x * x_res_ + x_min_);
	y = -1. * (y * y_res_ + y_min_);

	double alpha, delta;

	switch(proj_) {
	case Proj0: {
		delta = delta0_ - y;
		alpha = x / COS(delta / rad) + alpha0_;
		break;
	}
	case Proj9:
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
		delta = -ASIN(y * pv_[0]) * rad;
		alpha = x + alpha0_;
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
FlatSkyProjection::angle_to_xy(double alpha, double delta) const
{
	static const double circ = 360 * deg;
	static const double halfcirc = 180 * deg;
	double dalpha = alpha - alpha0_;

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
	case Proj9:
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
		y = -SIN(delta / rad) / pv_[0];
		break;
	}
	default:
		log_fatal("Proj %d not implemented", proj_);
		break;
	}

	x *= -1;
	y *= -1;

	x = (x - x_min_) / x_res_;
	y = (y - y_min_) / y_res_;

	return {x, y};
}

std::vector<double>
FlatSkyProjection::pixel_to_angle(long pixel, bool wrap_alpha) const
{
	if (pixel < 0 || pixel >= xpix_ * ypix_)
		return std::vector<double>(2, 0);
	std::vector<double> xy = pixel_to_xy(pixel);
	return xy_to_angle(xy[0], xy[1], wrap_alpha);
}

long
FlatSkyProjection::angle_to_pixel(double alpha, double delta) const
{
	std::vector<double> xy = angle_to_xy(alpha, delta);
	return xy_to_pixel(xy[0], xy[1]);
}

void FlatSkyProjection::get_rebin_angles(long pixel, size_t scale,
    std::vector<double> & alphas, std::vector<double> & deltas,
    bool wrap_alpha) const
{
	alphas = std::vector<double>(scale * scale);
	deltas = std::vector<double>(scale * scale);

	if (pixel < 0 || pixel >= xpix_ * ypix_) {
		log_debug("Point lies outside of pixel grid\n");
		return;
	}

	std::vector <double> xy = pixel_to_xy(pixel);
	double x0 = xy[0] - 0.5;
	double y0 = xy[1] - 0.5;

	size_t count = 0;
	for (size_t j = 0; j < scale; j++) {
		double y = y0 + (j + 0.5) / (double) scale;
		for (size_t i = 0; i < scale; i++) {
			double x = x0 + (i + 0.5) / (double) scale;
			std::vector<double> ang = xy_to_angle(x, y, wrap_alpha);
			alphas[count] = ang[0];
			deltas[count] = ang[1];
			count++;
		}
	}
}

void FlatSkyProjection::get_interp_pixels_weights(double alpha, double delta,
    std::vector<long> & pixels, std::vector<double> & weights) const
{
	std::vector<double> xy = angle_to_xy(alpha, delta);
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

FlatSkyProjection FlatSkyProjection::rebin(size_t scale) const
{
	FlatSkyProjection fp(*this);

	if (scale <= 1)
		return fp;

	fp.xpix_ /= scale;
	fp.ypix_ /= scale;
	fp.set_res(y_res_ * scale, x_res_ * scale);

	return fp;
}


G3_SPLIT_SERIALIZABLE_CODE(FlatSkyProjection);

PYBINDINGS("coordinateutils")
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
}
