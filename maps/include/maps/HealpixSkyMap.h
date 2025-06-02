#ifndef _MAPS_HEALPIXSKYMAP_H
#define _MAPS_HEALPIXSKYMAP_H

#include <G3Frame.h>
#include <G3Logging.h>

#include <vector>
#include <string>

#include <maps/G3SkyMap.h>
#include <maps/HealpixSkyMapInfo.h>

template <typename T> class SparseMapData;


class HealpixSkyMap : public G3FrameObject, public G3SkyMap {
public:
	// Construct a Healpix map with given nside, units, and coordinates.
	HealpixSkyMap(size_t nside,
	    bool weighted = true,
	    bool nested = false,
	    MapCoordReference coord_ref = MapCoordReference::Equatorial,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    G3SkyMap::MapPolType pol_type = G3SkyMap::None,
	    bool shift_ra = false,
	    G3SkyMap::MapPolConv pol_conv = G3SkyMap::ConvNone);

	HealpixSkyMap(const HealpixSkyMapInfo &info,
	    bool weighted = true,
	    MapCoordReference coord_ref = MapCoordReference::Equatorial,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    G3SkyMap::MapPolType pol_type = G3SkyMap::None,
	    G3SkyMap::MapPolConv pol_conv = G3SkyMap::ConvNone);

	HealpixSkyMap();
	HealpixSkyMap(const HealpixSkyMap& fm);

	~HealpixSkyMap();

	double at(size_t i) const override;
	double &operator [] (size_t i) override;
	double *data();

	// +
	virtual G3SkyMap &operator+=(const G3SkyMap &rhs) override;
	virtual G3SkyMap &operator+=(double rhs) override;

	// -
	virtual G3SkyMap &operator-=(const G3SkyMap &rhs) override;
	virtual G3SkyMap &operator-=(double rhs) override;

	// *
	virtual G3SkyMap &operator*=(const G3SkyMap &rhs) override;
	virtual G3SkyMap &operator*=(const G3SkyMapMask &rhs) override;
	virtual G3SkyMap &operator*=(double rhs) override;

	// /
	virtual G3SkyMap &operator/=(const G3SkyMap &rhs) override;
	virtual G3SkyMap &operator/=(double rhs) override;

	template <class A> void load(A &ar, unsigned v);
	template <class A> void save(A &ar, unsigned v) const;
	virtual G3SkyMapPtr Clone(bool copy_data = true) const override;
	std::string Description() const override;

	std::vector<size_t> shape() const override;
	size_t NpixAllocated() const override;
	size_t NpixNonZero() const override;
	bool IsCompatible(const G3SkyMap & other) const override;
	void NonZeroPixels(std::vector<uint64_t> &indices,
	    std::vector<double> &data) const; // Iterators better?
	void ApplyMask(const G3SkyMapMask &mask, bool inverse=false) override;

	size_t nside() const {return info_.nside();}
	bool nested() const {return info_.nested();}
	double res() const;

	size_t AngleToPixel(double alpha, double delta) const override;
	std::vector<double> PixelToAngle(size_t pixel) const override;
	size_t QuatToPixel(const Quat &q) const override;
	Quat PixelToQuat(size_t pixel) const override;

	G3VectorQuat GetRebinQuats(size_t pixel, size_t scale) const override;
	void GetInterpPixelsWeights(const Quat &q, std::vector<uint64_t> & pixels,
	    std::vector<double> & weights) const override;

	std::vector<uint64_t> QueryDisc(const Quat &q, double radius) const override;

	G3SkyMapPtr Rebin(size_t scale, bool norm = true) const override;

	void ConvertToDense() override;
	void ConvertToRingSparse();
	void ConvertToIndexedSparse();
	bool IsDense() const override { return (dense_ != NULL); }
	bool IsRingSparse() const { return (ring_sparse_ != NULL); }
	bool IsIndexedSparse() const { return (indexed_sparse_ != NULL); }
	void Compact(bool zero_nans = false) override;

	bool IsRaShifted() const { return info_.shifted(); }
	void SetShiftRa(bool shift);

	class const_iterator {
	public:
		typedef std::pair<uint64_t, double> value_type;
		typedef value_type & reference;
		typedef value_type * pointer;

		const_iterator(const HealpixSkyMap &map, bool begin);
		const_iterator(const const_iterator &iter);

		bool operator==(const const_iterator & other) const {
			return map_.ring_sparse_ ? (j_ == other.j_ && k_ == other.k_) :
			    index_ == other.index_;
		}
		bool operator!=(const const_iterator & other) const {
			return map_.ring_sparse_ ? (j_ != other.j_ || k_ != other.k_) :
			    index_ != other.index_;
		}

		reference operator*() { return value_; };
		pointer operator->() { return &value_; };

		const_iterator operator++();
		const_iterator operator++(int) { const_iterator i = *this; ++(*this); return i; }

	private:
		uint64_t index_;
		value_type value_;
		const HealpixSkyMap &map_;
		std::unordered_map<uint64_t, double>::const_iterator it_indexed_sparse_;
		std::vector<double>::const_iterator it_dense_;
		size_t j_, k_;

		void set_value();
	};

	const_iterator begin() const { return const_iterator(*this, true); };
	const_iterator end() const { return const_iterator(*this, false); };

private:
	HealpixSkyMapInfo info_;
	std::vector<double> *dense_;
	SparseMapData<double> *ring_sparse_;
	std::unordered_map<uint64_t, double> *indexed_sparse_;

	SET_LOGGER("HealpixSkyMap");
};

G3_POINTERS(HealpixSkyMap);
G3_SPLIT_SERIALIZABLE(HealpixSkyMap, 3);

#endif //_MAPS_HEALPIXSKYMAP_H

