#ifndef _COORDINATEUTILS_FLATSKYMAP_H
#define _COORDINATEUTILS_FLATSKYMAP_H

#include <G3Frame.h>
#include <G3Map.h>
#include <G3Logging.h>

#include <vector>
#include <string>

#include <coordinateutils/G3SkyMap.h>
#include <coordinateutils/flatskyprojection.h>

class FlatSkyMap : public G3SkyMap {
public:
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

	FlatSkyMap(const FlatSkyProjection & fp,
	    MapCoordReference coord_ref = MapCoordReference::Equatorial,
	    bool is_weighted = true,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    G3SkyMap::MapPolType pol_type = MapPolType::None);

	FlatSkyMap();
	FlatSkyMap(const FlatSkyMap & fm);

	template <class A> void load(A &ar, unsigned v);
	template <class A> void save(A &ar, unsigned v) const;
	virtual G3SkyMapPtr Clone(bool copy_data = true) const override;
	std::string Description() const override;

	bool IsCompatible(const G3SkyMap & other) const override;

	void set_proj(MapProjection proj);
	void set_alpha_center(double alpha);
	void set_delta_center(double delta);
	void set_center(double alpha, double delta);
	void set_xres(double res);
	void set_yres(double res);
	void set_res(double res);

	MapProjection proj() const;
	double alpha_center() const;
	double delta_center() const;
	double xres() const;
	double yres() const;
	double res() const;

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
	FlatSkyProjection proj_info; // projection parameters and functions

	SET_LOGGER("FlatSkyMap");
};

G3_POINTERS(FlatSkyMap);

namespace cereal {
  template <class A> struct specialize<A, FlatSkyMap, cereal::specialization::member_load_save> {};
}

G3_SERIALIZABLE(FlatSkyMap, 2);

#endif //_COORDINATEUTILS_FLATSKYMAP_H

