#ifndef _COORDINATEUTILS_FLATSKYMAP_H
#define _COORDINATEUTILS_FLATSKYMAP_H

#include <G3Frame.h>
#include <G3Logging.h>

#include <vector>
#include <string>

#include <coordinateutils/G3SkyMap.h>
#include <coordinateutils/flatskyprojection.h>

class DenseMapData;
class SparseMapData;

class FlatSkyMap : public G3FrameObject, public G3SkyMap {
public:
	// Construct an X by Y pixel flat map with pixel width res and the given
	// units, center, and coordinate system. If x_res is set to something
	// non-zero, will set the X resolution to a different number than res,
	// creating a map with rectangular pixels.
	FlatSkyMap(size_t x_len, size_t y_len, double res,
 	    bool is_weighted = true,
	    MapProjection proj = MapProjection::ProjNone, 
	    double alpha_center = 0, double delta_center = 0, 
	    MapCoordReference coord_ref = MapCoordReference::Equatorial,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    G3SkyMap::MapPolType pol_type = MapPolType::None,
	    double x_res = 0, /* if different from res */
	    double x_center = 0.0 / 0.0, double y_center = 0.0 / 0.0);

	// Constructor from a numpy array
	FlatSkyMap(boost::python::object v, double res, 
	    bool is_weighted = true,
	    MapProjection proj = MapProjection::ProjNone, 
	    double alpha_center = 0, double delta_center = 0, 
	    MapCoordReference coord_ref = MapCoordReference::Equatorial,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    G3SkyMap::MapPolType pol_type = MapPolType::None,
	    double x_res = 0, double x_center = 0.0 / 0.0,
	    double y_center = 0.0 / 0.0);

	FlatSkyMap(const FlatSkyProjection & fp,
	    MapCoordReference coord_ref = MapCoordReference::Equatorial,
	    bool is_weighted = true,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    G3SkyMap::MapPolType pol_type = MapPolType::None);

	FlatSkyMap();
	FlatSkyMap(const FlatSkyMap & fm);

	~FlatSkyMap();

        double operator [] (size_t i) const override;
        double &operator [] (size_t i) override;

        double operator () (size_t x, size_t y) const;
        double &operator () (size_t x, size_t y);
        double at(size_t x, size_t y) const { return (*this)(x, y); }

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
	    std::vector<double> &data) const;

	void set_proj(MapProjection proj);
	void set_alpha_center(double alpha);
	void set_delta_center(double delta);
	void set_x_center(double y);
	void set_y_center(double y);
	void set_xres(double res);
	void set_yres(double res);
	void set_res(double res);

	MapProjection proj() const;
	double alpha_center() const;
	double delta_center() const;
	double x_center() const;
	double y_center() const;
	double xres() const;
	double yres() const;
	double res() const;

	size_t angle_to_pixel(double alpha, double delta) const override;
	std::vector<double> pixel_to_angle(size_t pixel) const override;
        std::vector<double> pixel_to_angle(size_t x_pix, size_t y_pix) const;
	std::vector<double> pixel_to_angle_wrap_ra(size_t pixel) const;
	std::vector<double> angle_to_xy(double alpha, double delta) const;
	std::vector<double> xy_to_angle(double x, double y) const;

	void get_rebin_angles(long pixel, size_t scale,
	    std::vector<double> & alphas, std::vector<double> & deltas) const override;
	void get_interp_pixels_weights(double alpha, double delta,
	    std::vector<long> & pixels, std::vector<double> & weights) const override;

	G3SkyMapPtr Rebin(size_t scale, bool norm = true) const override;

	void ConvertToDense();
	void ConvertToSparse();
	bool IsDense() const { return (dense_ != NULL); }

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
	virtual void init_from_v1_data(std::vector<size_t>, const std::vector<double> &) override;

private:
	FlatSkyProjection proj_info; // projection parameters and functions

	DenseMapData *dense_;
	SparseMapData *sparse_;
	uint64_t xpix_, ypix_;

	SET_LOGGER("FlatSkyMap");
};

G3_POINTERS(FlatSkyMap);

namespace cereal {
  template <class A> struct specialize<A, FlatSkyMap, cereal::specialization::member_load_save> {};
}

G3_SERIALIZABLE(FlatSkyMap, 3);

#endif //_COORDINATEUTILS_FLATSKYMAP_H

