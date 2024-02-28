#ifndef _MAPS_FLATSKYMAP_H
#define _MAPS_FLATSKYMAP_H

#include <G3Frame.h>
#include <G3Logging.h>

#include <vector>
#include <string>

#include <maps/G3SkyMap.h>
#include <maps/FlatSkyProjection.h>

class DenseMapData;
template <typename T> class SparseMapData;

class FlatSkyMap : public G3FrameObject, public G3SkyMap {
public:
	// Construct an X by Y pixel flat map with pixel width res and the given
	// units, center, and coordinate system. If x_res is set to something
	// non-zero, will set the X resolution to a different number than res,
	// creating a map with rectangular pixels.
	FlatSkyMap(size_t x_len, size_t y_len, double res,
	    bool weighted = true,
	    MapProjection proj = MapProjection::ProjNone, 
	    double alpha_center = 0, double delta_center = 0, 
	    MapCoordReference coord_ref = MapCoordReference::Equatorial,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    G3SkyMap::MapPolType pol_type = G3SkyMap::None,
	    double x_res = 0, /* if different from res */
	    double x_center = 0.0 / 0.0, double y_center = 0.0 / 0.0,
	    bool flat_pol = false,
	    G3SkyMap::MapPolConv pol_conv = G3SkyMap::ConvNone);

	// Constructor from a numpy array
	FlatSkyMap(boost::python::object v, double res, 
	    bool weighted = true,
	    MapProjection proj = MapProjection::ProjNone, 
	    double alpha_center = 0, double delta_center = 0, 
	    MapCoordReference coord_ref = MapCoordReference::Equatorial,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    G3SkyMap::MapPolType pol_type = G3SkyMap::None,
	    double x_res = 0, double x_center = 0.0 / 0.0,
	    double y_center = 0.0 / 0.0, bool flat_pol = false,
	    G3SkyMap::MapPolConv pol_conv = G3SkyMap::ConvNone);

	FlatSkyMap(const FlatSkyProjection & fp,
	    MapCoordReference coord_ref = MapCoordReference::Equatorial,
	    bool weighted = true,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    G3SkyMap::MapPolType pol_type = G3SkyMap::None,
	    bool flat_pol = false,
	    G3SkyMap::MapPolConv pol_conv = G3SkyMap::ConvNone);

	FlatSkyMap();
	FlatSkyMap(const FlatSkyMap & fm);

	~FlatSkyMap();

	double at(size_t i) const override;
	double &operator [] (size_t i) override;

	double at(size_t x, size_t y) const;
	double operator () (size_t x, size_t y) const {	return this->at(x, y); };
	double &operator () (size_t x, size_t y);

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
	virtual void FillFromArray(boost::python::object v) override;
	virtual G3SkyMapPtr Clone(bool copy_data = true) const override;
	std::string Description() const override;

	std::vector<size_t> shape() const override;
	size_t NpixAllocated() const override;
	size_t NpixNonZero() const override;
	bool IsCompatible(const G3SkyMap & other) const override;
	void NonZeroPixels(std::vector<uint64_t> &indices,
	    std::vector<double> &data) const;
	void ApplyMask(const G3SkyMapMask &mask, bool inverse=false) override;

	void SetProj(MapProjection proj);
	void SetAlphaCenter(double alpha);
	void SetDeltaCenter(double delta);
	void SetXCenter(double y);
	void SetYCenter(double y);
	void SetXRes(double res);
	void SetYRes(double res);
	void SetRes(double res);

	MapProjection proj() const;
	double alpha_center() const;
	double delta_center() const;
	double x_center() const;
	double y_center() const;
	double xres() const;
	double yres() const;
	double res() const;

	size_t AngleToPixel(double alpha, double delta) const override;
	std::vector<double> PixelToAngle(size_t pixel) const override;
	std::vector<double> PixelToAngle(size_t x_pix, size_t y_pix) const;
	std::vector<double> AngleToXY(double alpha, double delta) const;
	std::vector<double> XYToAngle(double x, double y) const;
	size_t XYToPixel(double x, double y) const;
	std::vector<double> PixelToXY(size_t pixel) const;
	std::vector<double> QuatToXY(quat q) const;
	quat XYToQuat(double x, double y) const;
	size_t QuatToPixel(quat q) const override;
	quat PixelToQuat(size_t pixel) const override;

	std::vector<double> PixelToAngleGrad(size_t pixel, double h=0.001) const;

	G3VectorQuat GetRebinQuats(size_t pixel, size_t scale) const override;
	void GetInterpPixelsWeights(quat q, std::vector<size_t> & pixels,
	    std::vector<double> & weights) const override;

	std::vector<size_t> QueryDisc(quat q, double radius) const override;

	G3SkyMapPtr Rebin(size_t scale, bool norm = true) const override;
	G3SkyMapPtr ExtractPatch(size_t x0, size_t y0, size_t width, size_t height,
	    double fill = 0) const;
	void InsertPatch(const FlatSkyMap &patch, bool ignore_zeros = false);
	G3SkyMapPtr Reshape(size_t width, size_t height, double fill = 0) const;

	void ConvertToDense() override;
	void ConvertToSparse();
	bool IsDense() const override { return (dense_ != NULL); }
	void Compact(bool zero_nans = false) override;

	bool IsPolFlat() const { return flat_pol_; }
	void SetFlatPol(bool flat) { flat_pol_ = flat; }

	class const_iterator {
	public:
		typedef std::pair<uint64_t, double> value_type;
		typedef value_type & reference;
		typedef value_type * pointer;

		const_iterator(const FlatSkyMap &map, bool begin);

		bool operator==(const const_iterator & other) const {
			return ((x_ == other.x_) && (y_ == other.y_));
		}
		bool operator!=(const const_iterator & other) const {
			return ((x_ != other.x_) || (y_ != other.y_));;
		}

		reference operator*() { return value_; };
		pointer operator->() { return &value_; };

		const_iterator operator++();
		const_iterator operator++(int) { const_iterator i = *this; ++(*this); return i; }

	private:
		size_t x_, y_;
		value_type value_;
		const FlatSkyMap &map_;

		void set_value() {
			value_.first = x_ + y_ * map_.xpix_;
			value_.second = map_.at(x_, y_);
		}
	};

	const_iterator begin() const { return const_iterator(*this, true); };
	const_iterator end() const { return const_iterator(*this, false); };

protected:
	virtual void InitFromV1Data(std::vector<size_t>, const std::vector<double> &) override;

private:
	FlatSkyProjection proj_info; // projection parameters and functions

	DenseMapData *dense_;
	SparseMapData<double> *sparse_;
	uint64_t xpix_, ypix_;
	bool flat_pol_;

	SET_LOGGER("FlatSkyMap");
};

G3_POINTERS(FlatSkyMap);

namespace cereal {
	template <class A> struct specialize<A, FlatSkyMap, cereal::specialization::member_load_save> {};
}

G3_SERIALIZABLE(FlatSkyMap, 4);

#endif //_MAPS_FLATSKYMAP_H

