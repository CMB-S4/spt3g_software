#ifndef _MAPS_FLATSKYPROJECTION_H
#define _MAPS_FLATSKYPROJECTION_H

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
	ProjZenithalEquidistant = 3,
	ProjARC = 3,
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
			  double x_res = 0.0 / 0.0,
			  MapProjection proj = MapProjection::ProjNone,
			  double x_center = 0.0 / 0.0, double y_center = 0.0 / 0.0);

	FlatSkyProjection();
	FlatSkyProjection(const FlatSkyProjection & fp);

	void initialize(size_t xpix, size_t ypix, double res,
	    double alpha_center = 0, double delta_center = 0, double x_res = 0.0 / 0.0,
	    MapProjection proj = MapProjection::ProjNone,
	    double x_center = 0.0 / 0.0, double y_center = 0.0 / 0.0);

	template <class A> void load(A &ar, unsigned v);
	template <class A> void save(A &ar, unsigned v) const;
	bool IsCompatible(const FlatSkyProjection & other) const;
	std::string Description() const override;

	void SetProj(MapProjection proj);
	void SetAlphaCenter(double alpha);
	void SetDeltaCenter(double delta);
	void SetAngleCenter(double alpha, double delta);
	void SetXYCenter(double x, double y);
	void SetXCenter(double x);
	void SetYCenter(double y);
	void SetXRes(double res);
	void SetYRes(double res);
	void SetRes(double res, double x_res = 0.0 / 0.0);

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

	size_t XYToPixel(double x, double y) const;
	std::vector<double> PixelToXY(size_t pixel) const;
	std::vector<double> XYToAngle(double x, double y) const;
	std::vector<double> AngleToXY(double alpha, double delta) const;
	std::vector<double> PixelToAngle(size_t pixel) const;
	size_t AngleToPixel(double alpha, double delta) const;
	std::vector<double> QuatToXY(const Quat &q) const;
	Quat XYToQuat(double x, double y) const;
	size_t QuatToPixel(const Quat &q) const;
	Quat PixelToQuat(size_t pixel) const;

	std::vector<double> XYToAngleGrad(double x, double y, double h=0.001) const;
	std::vector<double> PixelToAngleGrad(size_t pixel, double h=0.001) const;

	size_t RebinPixel(size_t pixel, size_t scale) const;

	G3VectorQuat GetRebinQuats(size_t pixel, size_t scale) const;
	void GetInterpPixelsWeights(const Quat &q, std::vector<size_t> & pixels,
	    std::vector<double> & weights) const;

	std::vector<size_t> QueryDisc(const Quat &q, double radius) const;

	FlatSkyProjection Rebin(size_t scale, double x_center = 0.0 / 0.0,
	    double y_center = 0.0 / 0.0) const;
	FlatSkyProjection OverlayPatch(double x0, double y0, size_t width,
	    size_t height) const;
	std::vector<double> GetPatchCenter(const FlatSkyProjection &proj) const;

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
	bool cyl_;
	double sindelta0_;
	double cosdelta0_;
	Quat q0_;

	SET_LOGGER("FlatSkyProjection");
};

G3_POINTERS(FlatSkyProjection);
G3_SPLIT_SERIALIZABLE(FlatSkyProjection, 4);

#endif //#ifndef _MAPS_FLATSKYPROJECTION_H
