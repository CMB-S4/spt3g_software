#ifndef _COORDINATEUTILS_FLATSKYMAP_H
#define _COORDINATEUTILS_FLATSKYMAP_H

#include <G3Frame.h>
#include <G3Map.h>
#include <G3Logging.h>
#include <G3SkyMap.h>

#include <vector>
#include <string>

#include <coordinateutils/flatskyprojection.h>

class FlatSkyMap : public G3SkyMap {
public:
	MapProjection proj;
	double alpha_center;    // Map center position in X
	double delta_center;    // Map center in Y
	double res;             // Pixel height (and width if x_res is 0)
	double x_res;           // Pixel width if different from zero
	
	// Construct an X by Y pixel flat map with pixel width res and the given
	// units, center, and coordinate system. If x_res is set to something
	// non-zero, will set the X resolution to a different number than res,
	// creating a map with rectangular pixels.
	FlatSkyMap(int x_len, int y_len, double res, 
 	    bool is_weighted = true,
	    MapProjection proj = MapProjection::ProjNone, 
	    double alpha_center = 0, double delta_center = 0, 
	    MapCoordReference coord_ref = MapCoordReference::Equatorial,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    G3SkyMap::MapPolType pol_type = MapPolType::None,
	    double x_res = 0 /* if different from res */);

	// Constructor from a numpy array
	FlatSkyMap(boost::python::object v, double res, 
	    bool is_weighted = true,
	    MapProjection proj = MapProjection::ProjNone, 
	    double alpha_center = 0, double delta_center = 0, 
	    MapCoordReference coord_ref = MapCoordReference::Equatorial,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    G3SkyMap::MapPolType pol_type = MapPolType::None,
	    double x_res = 0);

	FlatSkyMap();
	FlatSkyMap(const FlatSkyMap & fm);

	template <class A> void serialize(A &ar, const unsigned u);
	virtual G3SkyMapPtr Clone(bool copy_data = true) const override;
	std::string Description() const override;

	bool IsCompatible(const G3SkyMap & other) const override;

	size_t angle_to_pixel(double alpha, double delta) const override;
	std::vector<double> pixel_to_angle(size_t pixel) const override;
	std::vector<double> pixel_to_angle_wrap_ra(size_t pixel) const;
	std::vector<double> angle_to_xy(double alpha, double delta) const;
	std::vector<double> xy_to_angle(double x, double y) const;

	void get_rebin_angles(long pixel, size_t scale,
	    std::vector<double> & alphas, std::vector<double> & deltas) const override;
	void get_interp_pixels_weights(double alpha, double delta,
	    std::vector<long> & pixels, std::vector<double> & weights) const override;

	G3SkyMapPtr rebin(size_t scale) const override;

private:
	FlatSkyProjection proj_info_;

	SET_LOGGER("FlatSkyMap");
};

G3_POINTERS(FlatSkyMap);
G3_SERIALIZABLE(FlatSkyMap, 1);

#endif //_COORDINATEUTILS_FLATSKYMAP_H

