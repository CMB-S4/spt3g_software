#include <pybindings.h>

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
    double alpha_center, double delta_center, double x_res, MapProjection proj) :
      alpha0_(alpha_center), delta0_(delta_center), xpix_(xpix), ypix_(ypix),
      x_res_(x_res > 0 ? x_res: res), y_res_(res), proj_(proj)
{
	x_min_ = -0.5 * (xpix_ - 1) * x_res_;
	y_min_ = -0.5 * (ypix_ - 1) * y_res_;
	sindelta0_ = SIN(delta0_ / rad);
	cosdelta0_ = COS(delta0_ / rad);
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
		x /= rad;
		double h = (delta0_ - y) / rad;
		double k = sqrt(1. - x * x);
		alpha = ATAN2(x, k * COS(h)) * rad + alpha0_;
		delta = ASIN(k * SIN(h)) * rad;
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
		dalpha /= rad;
		delta /= rad;
		double cosdelta = COS(delta);
		y = 90 * deg + delta0_ - ATAN2(COS(dalpha) * cosdelta, -SIN(delta)) * rad;
		x = cosdelta * SIN(dalpha) * rad;
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
