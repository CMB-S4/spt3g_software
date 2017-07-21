#ifndef _COORDINATEUTILS_FLATSKYMAP_H
#define _COORDINATEUTILS_FLATSKYMAP_H

#include <G3Frame.h>
#include <G3Map.h>
#include <G3Logging.h>
#include <G3SkyMap.h>

#include <vector>
#include <string>

// Defined map projections
enum MapProjection {
	// Standard projections
	ProjSansonFlamsteed = 0,
	ProjCAR = 1,
	ProjSIN = 2,
	ProjStereographic = 4,
	ProjLambertAzimuthalEqualArea = 5,
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
	    G3Timestream::TimestreamUnits u = G3Timestream::Kcmb,
	    G3SkyMap::MapPolType pol_type = MapPolType::None,
	    double x_res = 0 /* if different from res */);

	// Constructor from a numpy array
	FlatSkyMap(boost::python::object v, double res, 
	    bool is_weighted = true,
	    MapProjection proj = MapProjection::ProjNone, 
	    double alpha_center = 0, double delta_center = 0, 
	    MapCoordReference coord_ref = MapCoordReference::Equatorial,
	    G3Timestream::TimestreamUnits u = G3Timestream::Kcmb,
	    G3SkyMap::MapPolType pol_type = MapPolType::None,
	    double x_res = 0);

	FlatSkyMap();
	FlatSkyMap(const FlatSkyMap & fm);

	template <class A> void serialize(A &ar, const unsigned u);
	virtual G3SkyMapPtr Clone(bool copy_data) const override;
	std::string Description() const override;

	// Arithmetic operations:

	// +
	FlatSkyMap & operator+=(const FlatSkyMap & rhs);
	FlatSkyMap & operator+=(double rhs);
	FlatSkyMap operator+(const FlatSkyMap & rhs);
	FlatSkyMap operator+(double rhs);

	// -
	FlatSkyMap & operator-=(const FlatSkyMap & rhs);
	FlatSkyMap & operator-=(double rhs);
	FlatSkyMap operator-(const FlatSkyMap & rhs);
	FlatSkyMap operator-(double rhs);

	// *
	FlatSkyMap & operator*=(const FlatSkyMap & rhs);
	FlatSkyMap & operator*=(double rhs);
	FlatSkyMap operator*(const FlatSkyMap & rhs);
	FlatSkyMap operator*(double rhs);

	// /
	FlatSkyMap & operator/=(const FlatSkyMap & rhs);
	FlatSkyMap & operator/=(double rhs);
	FlatSkyMap operator/(const FlatSkyMap & rhs);
	FlatSkyMap operator/(double rhs);

	std::vector<int> angles_to_pixels(const std::vector<double> & alphas, 
	    const std::vector<double> & deltas) const override;
	std::vector<double> pixel_to_angle(size_t pixel) const override;
	std::vector<double> pixel_to_angle_wrap_ra(size_t pixel) const;

private:
	SET_LOGGER("FlatSkyMap");
};

G3_POINTERS(FlatSkyMap);
G3_SERIALIZABLE(FlatSkyMap, 1);

#endif //_COORDINATEUTILS_FLATSKYMAP_H

