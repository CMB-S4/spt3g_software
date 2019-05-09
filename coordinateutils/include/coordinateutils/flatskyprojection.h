#ifndef _COORDINATEUTILS_FLATSKYPROJECTION_H
#define _COORDINATEUTILS_FLATSKYPROJECTION_H

#include <vector>

#include <G3Logging.h>

typedef double angle_t;

// Defined map projections
enum MapProjection {
	// Standard projections
	ProjSansonFlamsteed = 0,
	ProjSFL = 0,
	ProjPlateCarree = 1,
	ProjCAR = 1,
	ProjOrthographic = 2,
	ProjSIN = 2,
	ProjStereographic = 4,
	ProjSTG = 4,
	ProjLambertAzimuthalEqualArea = 5,
	ProjZEA = 5,
	ProjGnomonic = 6,
	ProjTAN = 6,
	ProjCylindricalEqualArea = 7,
	ProjCEA = 7,
	ProjBICEP = 9,

	ProjNone = 42,

	// Compatibility aliases for SPTpol
	Proj0 = 0,
	Proj1 = 1,
	Proj2 = 2,
	Proj3 = 3,
	Proj4 = 4,
	Proj5 = 5,
	Proj6 = 6,
	Proj7 = 7,
	Proj8 = 8,
	Proj9 = 9,
};

class FlatSkyProjection {
public:
	FlatSkyProjection(size_t xpix, size_t ypix, double res,
			  double alpha_center = 0, double delta_center = 0,
			  double x_res = 0,
			  MapProjection proj = MapProjection::ProjNone);

	long xy_to_pixel(double x, double y) const;
	std::vector<double> pixel_to_xy(long pixel) const;
	std::vector<double> xy_to_angle(double x, double y, bool wrap_alpha=false) const;
	std::vector<double> angle_to_xy(double alpha, double delta) const;
	std::vector<double> pixel_to_angle(long pixel, bool wrap_alpha=false) const;
	long angle_to_pixel(double alpha, double delta) const;

        void get_rebin_angles(long pixel, size_t scale,
	    std::vector<double> & alphas, std::vector<double> & deltas,
	    bool wrap_alpha=false) const;
	void get_interp_pixels_weights(double alpha, double delta,
	    std::vector<long> & pixels, std::vector<double> & weights) const;

private:
	// projection parameters
	double alpha0_;
	double delta0_;
	size_t xpix_;
	size_t ypix_;
	double x_res_;
	double y_res_;
	MapProjection proj_;

	// derived values
	double x_min_;
	double y_min_;
	double sindelta0_;
	double cosdelta0_;

	SET_LOGGER("FlatSkyProjection");
};

G3_POINTERS(FlatSkyProjection);

#endif //#ifndef _COORDINATEUTILS_FLATSKYPROJECTION_H
