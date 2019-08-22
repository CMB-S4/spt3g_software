#ifndef _COORDINATEUTILS_HEALPIXSKYMAP_H
#define _COORDINATEUTILS_HEALPIXSKYMAP_H

#include <G3Frame.h>
#include <G3Logging.h>

#include <vector>
#include <string>

#include <coordinateutils/G3SkyMap.h>
#include <coordinateutils/chealpix.h>

class SparseMapData;

class HealpixSkyMap : public G3FrameObject, public G3SkyMap {
public:
	// Construct an X by Y pixel flat map with pixel width res and the given
	// units, center, and coordinate system. If x_res is set to something
	// non-zero, will set the X resolution to a different number than res,
	// creating a map with rectangular pixels.
	HealpixSkyMap(size_t nside,
 	    bool is_weighted = true,
	    MapCoordReference coord_ref = MapCoordReference::Equatorial,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    G3SkyMap::MapPolType pol_type = MapPolType::None);

	// Constructor from a (dense) numpy array
	HealpixSkyMap(boost::python::object v,
	    bool is_weighted = true,
	    MapCoordReference coord_ref = MapCoordReference::Equatorial,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    G3SkyMap::MapPolType pol_type = MapPolType::None);

	// Constructor from a (sparse) numpy array
	HealpixSkyMap(boost::python::object indices, boost::python::object v,
	    size_t nside, bool is_weighted = true,
	    MapCoordReference coord_ref = MapCoordReference::Equatorial,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    G3SkyMap::MapPolType pol_type = MapPolType::None);

	HealpixSkyMap();
	HealpixSkyMap(const HealpixSkyMap& fm);

	~HealpixSkyMap();

        double operator [] (int i) const override;
        double &operator [] (int i) override;

	// XXX Usefully implement in-place operators for sparse maps

	template <class A> void load(A &ar, unsigned v);
	template <class A> void save(A &ar, unsigned v) const;
	virtual G3SkyMapPtr Clone(bool copy_data = true) const override;
	std::string Description() const override;

	std::vector<size_t> shape() const override;
	size_t npix_allocated() const override;
	bool IsCompatible(const G3SkyMap & other) const override;

	size_t nside() const;

	size_t angle_to_pixel(double alpha, double delta) const override;
	std::vector<double> pixel_to_angle(size_t pixel) const override;

	void get_rebin_angles(long pixel, size_t scale,
	    std::vector<double> & alphas, std::vector<double> & deltas) const override;
	void get_interp_pixels_weights(double alpha, double delta,
	    std::vector<long> & pixels, std::vector<double> & weights) const override;

	G3SkyMapPtr rebin(size_t scale) const override;

	void ConvertToDense();
	void ConvertToSparse();
	bool IsDense() const { return (dense_ != NULL); }

private:
	uint32_t nside_;
	bool is_nested_;
	std::vector<double> *dense_;
	SparseMapData *ring_sparse_;
	std::unordered_map<uint32_t, double> *indexed_sparse_;
	map_info *ring_info_;

	SET_LOGGER("HealpixSkyMap");
};

G3_POINTERS(HealpixSkyMap);

namespace cereal {
  template <class A> struct specialize<A, HealpixSkyMap, cereal::specialization::member_load_save> {};
}

G3_SERIALIZABLE(HealpixSkyMap, 1);

#endif //_COORDINATEUTILS_HEALPIXSKYMAP_H

