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
	// Construct a Healpix map with given nside, units, and coordinates.
	HealpixSkyMap(size_t nside,
 	    bool is_weighted = true,
 	    bool is_nested = false,
	    MapCoordReference coord_ref = MapCoordReference::Equatorial,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    G3SkyMap::MapPolType pol_type = MapPolType::None);

	// Constructor from a numpy array
	HealpixSkyMap(boost::python::object v,
	    bool is_weighted = true,
 	    bool is_nested = false,
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
	void NonZeroPixels(std::vector<uint32_t> &indices,
	    std::vector<double> &data) const; // Iterators better?

	size_t nside() const {return nside_;}
	bool nested() const {return is_nested_;}

	size_t angle_to_pixel(double alpha, double delta) const override;
	std::vector<double> pixel_to_angle(size_t pixel) const override;

	void get_rebin_angles(long pixel, size_t scale,
	    std::vector<double> & alphas, std::vector<double> & deltas) const override;
	void get_interp_pixels_weights(double alpha, double delta,
	    std::vector<long> & pixels, std::vector<double> & weights) const override;

	G3SkyMapPtr rebin(size_t scale) const override;

	void ConvertToDense();
	void ConvertToRingSparse();
	void ConvertToIndexedSparse();
	bool IsDense() const { return (dense_ != NULL); }
	bool IsRingSparse() const { return (indexed_sparse_ != NULL); }
	bool IsIndexedSparse() const { return (indexed_sparse_ != NULL); }

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

