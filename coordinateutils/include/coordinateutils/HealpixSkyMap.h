#ifndef _COORDINATEUTILS_HEALPIXSKYMAP_H
#define _COORDINATEUTILS_HEALPIXSKYMAP_H

#include <G3Frame.h>
#include <G3Logging.h>

#include <vector>
#include <string>

#include <coordinateutils/G3SkyMap.h>
#include <coordinateutils/chealpix.h>

class SparseMapData;
class SparseMapIterator;


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

        double operator [] (size_t i) const override;
        double &operator [] (size_t i) override;

	// +
	virtual G3SkyMap &operator+=(const G3SkyMap &rhs) override;
	virtual G3SkyMap &operator+=(double rhs) override;

	// -
	virtual G3SkyMap &operator-=(const G3SkyMap &rhs) override;
	virtual G3SkyMap &operator-=(double rhs) override;

	// *
	virtual G3SkyMap &operator*=(const G3SkyMap &rhs) override;
	virtual G3SkyMap &operator*=(double rhs) override;

	// /
	virtual G3SkyMap &operator/=(const G3SkyMap &rhs) override;
	virtual G3SkyMap &operator/=(double rhs) override;

	template <class A> void load(A &ar, unsigned v);
	template <class A> void save(A &ar, unsigned v) const;
	virtual G3SkyMapPtr Clone(bool copy_data = true) const override;
	std::string Description() const override;

	std::vector<size_t> shape() const override;
	size_t npix_allocated() const override;
	bool IsCompatible(const G3SkyMap & other) const override;
	void NonZeroPixels(std::vector<uint64_t> &indices,
	    std::vector<double> &data) const; // Iterators better?

	size_t nside() const {return nside_;}
	bool nested() const {return is_nested_;}

	size_t angle_to_pixel(double alpha, double delta) const override;
	std::vector<double> pixel_to_angle(size_t pixel) const override;

	void get_rebin_angles(long pixel, size_t scale,
	    std::vector<double> & alphas, std::vector<double> & deltas) const override;
	void get_interp_pixels_weights(double alpha, double delta,
	    std::vector<long> & pixels, std::vector<double> & weights) const override;

	G3SkyMapPtr Rebin(size_t scale, bool norm = true) const override;

	void ConvertToDense();
	void ConvertToRingSparse();
	void ConvertToIndexedSparse();
	bool IsDense() const { return (dense_ != NULL); }
	bool IsRingSparse() const { return (ring_sparse_ != NULL); }
	bool IsIndexedSparse() const { return (indexed_sparse_ != NULL); }

	class iterator {
	public:
		typedef std::pair<uint64_t, double> value_type;
		typedef value_type & reference;
		typedef value_type * pointer;

		iterator(const HealpixSkyMap &map, bool begin);
		iterator(const iterator &iter);

		bool operator==(const iterator & other) const {
			return (index_ == other.index_);
		}
		bool operator!=(const iterator & other) const {
			return (index_ != other.index_);
		}

		reference operator*() { return value_; };
		pointer operator->() { return &value_; };

		iterator operator++();
		iterator operator++(int) { iterator i = *this; ++(*this); return i; }

	private:
		uint64_t index_;
		value_type value_;
		const HealpixSkyMap &map_;
		std::unordered_map<uint64_t, double>::iterator it_indexed_sparse_;
		std::vector<double>::iterator it_dense_;
		size_t j_, k_;

		void set_value();
	};

	iterator begin() const { return iterator(*this, true); };
	iterator end() const { return iterator(*this, false); };

private:
	uint32_t nside_;
	size_t npix_;
	bool is_nested_;
	std::vector<double> *dense_;
	SparseMapData *ring_sparse_;
	std::unordered_map<uint64_t, double> *indexed_sparse_;
	map_info *ring_info_;

	SET_LOGGER("HealpixSkyMap");
};

G3_POINTERS(HealpixSkyMap);

namespace cereal {
  template <class A> struct specialize<A, HealpixSkyMap, cereal::specialization::member_load_save> {};
}

G3_SERIALIZABLE(HealpixSkyMap, 1);

#endif //_COORDINATEUTILS_HEALPIXSKYMAP_H

