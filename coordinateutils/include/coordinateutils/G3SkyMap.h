#ifndef _G3_SKYMAP_H
#define _G3_SKYMAP_H

#include <G3Frame.h>
#include <G3Timestream.h>
#include <vector>
#include <map>

#include <boost/python.hpp>

enum MapCoordReference {
	Local = 0,
	Equatorial = 1,
	Galactic = 2
};

/*
 * The following implements a generic container for a sky map. In general,
 * you will use subclasses of this object (see: FlatSkyMap and HealpixSkymap)
 * rather than this class directly. This class is meant to be a generic
 * interface for the map maker.
 */

class G3SkyMap {
public:
	// Following numerical values are important
	enum MapPolType {
		T = 0,
		Q = 1,
		U = 2,
		None = 7
	};
	
	G3SkyMap(MapCoordReference coords, size_t xpix, size_t ypix = 1,
	    bool isweighted = true,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    MapPolType pol_type = None) :
	    coord_ref(coords), units(u), pol_type(pol_type),
	    is_weighted(isweighted), xpix_(xpix), ypix_(ypix) {}
	virtual ~G3SkyMap() {};

	// Reimplement the following in subclasses
	virtual boost::shared_ptr<G3SkyMap> Clone(bool copy_data = true) const = 0;

	MapCoordReference coord_ref;
	G3Timestream::TimestreamUnits units;
	MapPolType pol_type;
	bool is_weighted;
	double overflow;

	void SetOverflow(double val) __attribute__((deprecated)) { 
		overflow = val; 
	}
	double GetOverflow(void) const __attribute__((deprecated)) {
		return overflow;
	}

	// Return a (modifiable) pixel value
	virtual double &operator [] (int i) = 0;
	virtual double operator [] (int i) const = 0;

	size_t xdim(void) const {
		return xpix_;
	}
	size_t ydim(void) const {
		return ypix_;
	}
	size_t size(void) const {
		return ypix_*xpix_;
	}

	virtual bool IsCompatible(const G3SkyMap & other) const {
		return ((xpix_ == other.xpix_) &&
			(ypix_ == other.ypix_) &&
			(coord_ref == other.coord_ref));
	}

	// Arithmetic operations:

#if 0
	// +
	G3SkyMap & operator+=(const G3SkyMap & rhs);
	G3SkyMap & operator+=(double rhs);
	G3SkyMap operator+(const G3SkyMap & rhs);
	G3SkyMap operator+(double rhs);

	// -
	G3SkyMap & operator-=(const G3SkyMap & rhs);
	G3SkyMap & operator-=(double rhs);
	G3SkyMap operator-(const G3SkyMap & rhs);
	G3SkyMap operator-(double rhs);

	// *
	G3SkyMap & operator*=(const G3SkyMap & rhs);
	G3SkyMap & operator*=(double rhs);
	G3SkyMap operator*(const G3SkyMap & rhs);
	G3SkyMap operator*(double rhs);

	// /
	G3SkyMap & operator/=(const G3SkyMap & rhs);
	G3SkyMap & operator/=(double rhs);
	G3SkyMap operator/(const G3SkyMap & rhs);
	G3SkyMap operator/(double rhs);
#endif

	// Pointing information
	virtual std::vector<int> angles_to_pixels(const std::vector<double> & alphas, 
	    const std::vector<double> & deltas) const;
	virtual void pixels_to_angles(const std::vector<int> & pixels,
	    std::vector<double> & alphas, std::vector<double> & deltas) const;

	virtual std::vector<double> pixel_to_angle(size_t pixel) const = 0;
	virtual std::vector<double> pixel_to_angle(size_t x_pix, size_t y_pix) const = 0;
	virtual size_t angle_to_pixel(double alpha, double delta) const  = 0;

	// Rebinning and interpolation
	virtual void get_rebin_angles(long pixel, size_t scale,
	    std::vector<double> & alphas, std::vector<double> & deltas) const = 0;
	virtual void get_interp_pixels_weights(double alpha, double delta,
	    std::vector<long> & pixels, std::vector<double> & weights) const = 0;
	double get_interp_precalc(const std::vector<long> & pixels,
	    const std::vector<double> & weights) const;
	double get_interp_value(double alpha, double delta) const;
	std::vector<double> get_interp_values(const std::vector<double> & alphas,
	    const std::vector<double> & deltas) const ;

	virtual boost::shared_ptr<G3SkyMap> rebin(size_t scale) const = 0;

protected:
	uint32_t xpix_, ypix_;
	virtual void init_from_v1_data(const std::vector<double> &) {
		throw std::runtime_error("Initializing from V1 not implemented");
	}

private:
	G3SkyMap() {} // Fake out for serialization
	template <class A> void serialize(A &ar, const unsigned v);
	friend class cereal::access;

	SET_LOGGER("G3SkyMap");
};

G3_POINTERS(G3SkyMap);

/*
 * This class implements a reference to another backing store of a 3x3
 * covariance matrix and is used to provide a matrix-based interface to
 * weight maps for the map maker.
 */
class PolCovMatrix {
public:
	PolCovMatrix(double &tt_, double &tq_, double &tu_, double &qq_,
	    double &qu_, double &uu_) : tt(tt_), tq(tq_), tu(tu_), qq(qq_),
	    qu(qu_), uu(uu_) {}
	PolCovMatrix(const PolCovMatrix &r) : tt(r.tt), tq(r.tq), tu(r.tu),
	    qq(r.qq), qu(r.qu), uu(r.uu) {}

	// Note: the default constructor uses an internal buffer. This lets
	// you initialize one of these statically and then use it in arithmetic.
	PolCovMatrix() : tt(backing[0]), tq(backing[1]),
	    tu(backing[2]), qq(backing[3]), qu(backing[4]), uu(backing[5]) {}

	double &tt, &tq, &tu, &qq, &qu, &uu;

	PolCovMatrix &operator +=(const PolCovMatrix &r) {
		tt += r.tt; tq += r.tq; tu += r.tu;
		qq += r.qq; qu += r.qu; uu += r.uu;
		return *this;
	}
	PolCovMatrix &operator -=(const PolCovMatrix &r) {
		tt -= r.tt; tq -= r.tq; tu -= r.tu;
		qq -= r.qq; qu -= r.qu; uu -= r.uu;
		return *this;
	}
	PolCovMatrix &operator =(const PolCovMatrix &r) {
		tt = r.tt; tq = r.tq; tu = r.tu;
		qq = r.qq; qu = r.qu; uu = r.uu;
		return *this;
	}

	PolCovMatrix &operator *=(const double r) {
		tt *= r; tq *= r; tu *= r;
		qq *= r; qu *= r; uu *= r;
		return *this;
	}

private:
	double backing[6];
};

class G3SkyMapWeights : public G3FrameObject {
public:
	enum WeightType {
		Wpol = 3,
		Wunpol = 4,
		None = 5
	};

	G3SkyMapWeights() : weight_type(None) {}

	// Instantiate weight maps based on the metadata of a reference map
	G3SkyMapWeights(G3SkyMapConstPtr ref_map, WeightType wt);

	G3SkyMapWeights(const G3SkyMapWeights &r);
	G3SkyMapPtr TT, TQ, TU, QQ, QU, UU;

	WeightType weight_type;

	PolCovMatrix operator [] (int i) {
		return PolCovMatrix((*TT)[i], (*TQ)[i], (*TU)[i], (*QQ)[i],
		    (*QU)[i], (*UU)[i]);
	}

	const PolCovMatrix at (int i) const {
		return PolCovMatrix((*TT)[i], (*TQ)[i], (*TU)[i], (*QQ)[i],
		    (*QU)[i], (*UU)[i]);
	}
	
#if 0
	void EnsureAllocated(void) {
		if (weight_type == Wunpol) {
			TT->EnsureAllocated();
		} else {
			TT->EnsureAllocated();
			TQ->EnsureAllocated();
			TU->EnsureAllocated();
			QQ->EnsureAllocated();
			QU->EnsureAllocated();
			UU->EnsureAllocated();
		}
	}
#endif

	boost::shared_ptr<G3SkyMapWeights> Clone(bool copy_data) const {
		if (copy_data)
			return boost::make_shared<G3SkyMapWeights>(*this);
		else
			return boost::make_shared<G3SkyMapWeights>(this->TT, this->weight_type);
	}
private:
	template <class A> void serialize(A &ar, const unsigned v);
	friend class cereal::access;
};

G3_POINTERS(G3SkyMapWeights);

G3_SERIALIZABLE(G3SkyMap, 2);
G3_SERIALIZABLE(G3SkyMapWeights, 2);

#endif

