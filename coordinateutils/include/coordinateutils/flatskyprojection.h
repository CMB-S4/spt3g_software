#ifndef _COORDINATEUTILS_FLATSKYPROJECTION_H
#define _COORDINATEUTILS_FLATSKYPROJECTION_H

#include <vector>

#include <G3Frame.h>
#include <G3Logging.h>

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

class FlatSkyProjection : public G3FrameObject {
public:
	FlatSkyProjection(size_t xpix, size_t ypix, double res,
			  double alpha_center = 0, double delta_center = 0,
			  double x_res = 0,
			  MapProjection proj = MapProjection::ProjNone,
			  double x_center = 0.0 / 0.0, double y_center = 0.0 / 0.0);

	FlatSkyProjection();
	FlatSkyProjection(const FlatSkyProjection & fp);

	void initialize(size_t xpix, size_t ypix, double res,
	    double alpha_center = 0, double delta_center = 0, double x_res = 0,
	    MapProjection proj = MapProjection::ProjNone,
	    double x_center = 0.0 / 0.0, double y_center = 0.0 / 0.0);

	template <class A> void load(A &ar, unsigned v);
	template <class A> void save(A &ar, unsigned v) const;
	bool IsCompatible(const FlatSkyProjection & other) const;
	std::string Description() const override;

	void set_proj(MapProjection proj);
	void set_alpha_center(double alpha);
	void set_delta_center(double delta);
	void set_angle_center(double alpha, double delta);
	void set_xy_center(double x, double y);
	void set_x_center(double x);
	void set_y_center(double y);
	void set_xres(double res);
	void set_yres(double res);
	void set_res(double res, double x_res=0);

	size_t xdim() const { return xpix_; };
	size_t ydim() const { return ypix_; };
	MapProjection proj() const { return proj_; };
	double alpha_center() const { return alpha0_; };
	double delta_center() const { return delta0_; };
	double x_center() const { return x0_; };
	double y_center() const { return y0_; };
	double xres() const { return x_res_; };
	double yres() const { return y_res_; };
	double res() const { return y_res_; };

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

	FlatSkyProjection rebin(size_t scale) const;

private:
	// projection parameters
	size_t xpix_;
	size_t ypix_;
	MapProjection proj_;
	double alpha0_;
	double delta0_;
	double x0_;
	double y0_;
	double x_res_;
	double y_res_;

	// derived values
	double x_min_;
	double y_min_;
	double sindelta0_;
	double cosdelta0_;
	std::vector<double> pv_;

	SET_LOGGER("FlatSkyProjection");
};

G3_POINTERS(FlatSkyProjection);

namespace cereal {
  template <class A> struct specialize<A, FlatSkyProjection, cereal::specialization::member_load_save> {};
}

G3_SERIALIZABLE(FlatSkyProjection, 3);

#endif //#ifndef _COORDINATEUTILS_FLATSKYPROJECTION_H
