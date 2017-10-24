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

class G3SkyMap : public G3FrameObject {
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
	    G3Timestream::TimestreamUnits u = G3Timestream::Kcmb,
	    MapPolType pol_type = None, 
	    bool pre_allocate_map = false) :
	    coord_ref(coords), units(u), pol_type(pol_type),
	    is_weighted(isweighted), data_(pre_allocate_map ? xpix*ypix+1 : 0),
	    xpix_(xpix), ypix_(ypix) {}

	// Instantiate from a numpy array
	G3SkyMap(boost::python::object v, MapCoordReference coords,
	    bool is_weighted, G3Timestream::TimestreamUnits units,
	    G3SkyMap::MapPolType pol_type);

	// Reimplement the following in subclasses
	virtual boost::shared_ptr<G3SkyMap> Clone(bool copy_data = true) const
	{
		if (copy_data)
			return boost::make_shared<G3SkyMap>(*this);
		else
			return boost::make_shared<G3SkyMap>(coord_ref, xpix_,
			    ypix_, is_weighted, units, pol_type, false);
	}

	MapCoordReference coord_ref;
	G3Timestream::TimestreamUnits units;
	MapPolType pol_type;
	bool is_weighted;

	void SetOverflow(double val) { 
		EnsureAllocated();
		*(data_.end() - 1) = val; 
	}
	double GetOverflow(void) const { 
		if (!IsAllocated()) return 0;
		return *(data_.end() - 1); 
	}

	// Return a (modifiable) pixel value
	double &operator [] (int i) {
		return data_[i];
	}
	double operator [] (int i) const {
		return data_[i];
	}
	const double &at (int i) const {
		return data_.at(i);
	}

	int pixat(int x, int y) const {
		return y*xpix_ + x;
	}

	size_t size(void) const {
		return (data_.size() == 0) ? 0 : (data_.size() - 1);
	}

	size_t xdim(void) const {
		return xpix_;
	}
	size_t ydim(void) const {
		return ypix_;
	}

	void EnsureAllocated(void) {
		if (!IsAllocated()) 
			data_.resize( xpix_ * ypix_ + 1, 0);
	}

	bool IsAllocated(void) const {
		return data_.size() > 0;
	}

	//pointing information
	virtual std::vector<int> angles_to_pixels(const std::vector<double> & alphas, 
					     const std::vector<double> & deltas) const { 
		return std::vector<int>();
	}

	virtual std::vector<double> pixel_to_angle(size_t pixel) const {
		return std::vector<double>();
	}

	std::vector<double> pixel_to_angle(size_t x_pix, size_t y_pix) const;
	size_t angle_to_pixel(double alpha, double delta) const;

protected:
	// Last bin is an overflow bin
	std::vector<double> data_;
	uint32_t xpix_, ypix_;

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

G3_SERIALIZABLE(G3SkyMap, 1);
G3_SERIALIZABLE(G3SkyMapWeights, 2);

#endif

