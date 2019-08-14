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

	G3SkyMap(MapCoordReference coords, bool isweighted = true,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    MapPolType pol_type = None) :
	    coord_ref(coords), units(u), pol_type(pol_type),
	    is_weighted(isweighted) {}
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

	virtual size_t size(void) const;  // total number of pixels
	virtual size_t npix_allocated(void) const = 0;  // stored in RAM
	virtual std::vector<size_t> shape(void) const = 0;  // map shape

	virtual bool IsCompatible(const G3SkyMap & other) const {
		if (shape().size() != other.shape().size())
			return false;
		for (size_t i = 0; i < shape().size(); i++)
			if (shape()[i] != other.shape()[i])
				return false;
		return (coord_ref == other.coord_ref);
	}

	// Arithmetic operations:

	// +
	virtual G3SkyMap &operator+=(const G3SkyMap &rhs);
	virtual G3SkyMap &operator+=(double rhs);

	// -
	virtual G3SkyMap &operator-=(const G3SkyMap &rhs);
	virtual G3SkyMap &operator-=(double rhs);

	// *
	virtual G3SkyMap &operator*=(const G3SkyMap &rhs);
	virtual G3SkyMap &operator*=(double rhs);

	// /
	virtual G3SkyMap &operator/=(const G3SkyMap &rhs);
	virtual G3SkyMap &operator/=(double rhs);

	// Pointing information
	virtual std::vector<int> angles_to_pixels(const std::vector<double> & alphas,
	    const std::vector<double> & deltas) const;
	virtual void pixels_to_angles(const std::vector<int> & pixels,
	    std::vector<double> & alphas, std::vector<double> & deltas) const;

	virtual std::vector<double> pixel_to_angle(size_t pixel) const = 0;
	virtual size_t angle_to_pixel(double alpha, double delta) const  = 0;

	// Rebinning and interpolation
	virtual void get_rebin_angles(long pixel, size_t scale,
	    std::vector<double> & alphas, std::vector<double> & deltas) const = 0;
	virtual void get_interp_pixels_weights(double alpha, double delta,
	    std::vector<long> & pixels, std::vector<double> & weights) const = 0;
	virtual double get_interp_precalc(const std::vector<long> & pixels,
	    const std::vector<double> & weights) const = 0;
	virtual double get_interp_value(double alpha, double delta) const;
	virtual std::vector<double> get_interp_values(const std::vector<double> & alphas,
	    const std::vector<double> & deltas) const;

	virtual boost::shared_ptr<G3SkyMap> rebin(size_t scale) const = 0;

protected:
	virtual void init_from_v1_data(const std::vector<size_t> &, const std::vector<double> &) {
		throw std::runtime_error("Initializing from V1 not implemented");
	}

private:
	G3SkyMap() {} // Fake out for serialization
	template <class A> void serialize(A &ar, const unsigned v);
	friend class cereal::access;

	SET_LOGGER("G3SkyMap");
};

G3_POINTERS(G3SkyMap);

class StokesVector;
class MuellerMatrix;

/*
 * This class implements a reference to another backing store of a 1x3
 * vector and is used to provide a matrix-based interface to
 * Stokes terms for the mapmaker.
 */
class StokesVector {
public:
	StokesVector(double &t_, double &q_, double &u_) :
	    t(t_), q(q_), u(u_) {}
	StokesVector(const StokesVector &r) : t(r.t), q(r.q), u(r.u) {}

	// Note: the default constructor uses an internal buffer. This lets
	// you initialize one of these statically and then use it in arithmetic.
	StokesVector() : t(backing[0]), q(backing[1]), u(backing[2]) {}

	double &t, &q, &u;

	StokesVector &operator +=(const StokesVector &r) {
		t += r.t; q += r.q; u += r.u;
		return *this;
	}
	StokesVector &operator -=(const StokesVector &r) {
		t -= r.t; q -= r.q; u -= r.u;
		return *this;
	}
	StokesVector &operator =(const StokesVector &r) {
		t = r.t; q = r.q; u = r.u;
		return *this;
	}

	StokesVector &operator *=(const double r) {
		t *= r; q *= r; u *= r;
		return *this;
	}

	StokesVector &operator /=(const MuellerMatrix &r);
	StokesVector operator /(const MuellerMatrix &r) const;

private:
	double backing[3];
};

/*
 * This class implements a reference to another backing store of a 3x3
 * covariance matrix and is used to provide a matrix-based interface to
 * weight maps for the map maker.
 */
class MuellerMatrix {
public:
	MuellerMatrix(double &tt_, double &tq_, double &tu_, double &qq_,
	    double &qu_, double &uu_) : tt(tt_), tq(tq_), tu(tu_), qq(qq_),
	    qu(qu_), uu(uu_) {}
	MuellerMatrix(const MuellerMatrix &r) : tt(r.tt), tq(r.tq), tu(r.tu),
	    qq(r.qq), qu(r.qu), uu(r.uu) {}

	// Note: the default constructor uses an internal buffer. This lets
	// you initialize one of these statically and then use it in arithmetic.
	MuellerMatrix() : tt(backing[0]), tq(backing[1]),
	    tu(backing[2]), qq(backing[3]), qu(backing[4]), uu(backing[5]) {}

	double &tt, &tq, &tu, &qq, &qu, &uu;

	MuellerMatrix &operator +=(const MuellerMatrix &r) {
		tt += r.tt; tq += r.tq; tu += r.tu;
		qq += r.qq; qu += r.qu; uu += r.uu;
		return *this;
	}
	MuellerMatrix &operator -=(const MuellerMatrix &r) {
		tt -= r.tt; tq -= r.tq; tu -= r.tu;
		qq -= r.qq; qu -= r.qu; uu -= r.uu;
		return *this;
	}
	MuellerMatrix &operator =(const MuellerMatrix &r) {
		tt = r.tt; tq = r.tq; tu = r.tu;
		qq = r.qq; qu = r.qu; uu = r.uu;
		return *this;
	}

	MuellerMatrix &operator *=(const double r) {
		tt *= r; tq *= r; tu *= r;
		qq *= r; qu *= r; uu *= r;
		return *this;
	}

	StokesVector operator *(const StokesVector &r) const {
		StokesVector s;
		s.t = tt * r.t + tq * r.q + tu * r.u;
		s.q = tq * r.t + qq * r.q + qu * r.u;
		s.u = tu * r.t + qu * r.q + uu * r.u;
		return s;
	}

	double det() const {
		return (tt * (qq * uu - qu * qu) -
			tq * (tq * uu - qu * tu) +
			tu * (tq * qu - qq * tu));
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

	MuellerMatrix operator [] (int i) {
		return MuellerMatrix((*TT)[i], (*TQ)[i], (*TU)[i], (*QQ)[i],
		    (*QU)[i], (*UU)[i]);
	}

	const MuellerMatrix at (int i) const {
		return MuellerMatrix((*TT)[i], (*TQ)[i], (*TU)[i], (*QQ)[i],
		    (*QU)[i], (*UU)[i]);
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

class G3SkyMapWithWeights : public G3FrameObject {
public:
	G3SkyMapWithWeights() {}

	G3SkyMapWithWeights(G3SkyMapConstPtr ref_map, bool isweighted = true,
	    bool ispolarized = true, std::string map_id = "");

	G3SkyMapWithWeights(G3SkyMapPtr T, G3SkyMapWeightsPtr weights,
	    bool isweighted = true, std::string map_id = "");
	G3SkyMapWithWeights(G3SkyMapPtr T, G3SkyMapPtr Q, G3SkyMapPtr U,
	    G3SkyMapWeightsPtr weights, bool isweighted = true,
	    std::string map_id = "");

	G3SkyMapWithWeights(const G3SkyMapWithWeights &r);

	std::string map_id;

	G3SkyMapPtr T, Q, U;
	G3SkyMapWeightsPtr weights;

	StokesVector operator [] (int i) {
		return StokesVector((*T)[i], (*Q)[i], (*U)[i]);
	}

	const StokesVector at (int i) const {
		return StokesVector((*T)[i], (*Q)[i], (*U)[i]);
	}

	bool IsWeighted() const {
		return !!weights;
	}

	bool IsPolarized() const {
		return !!Q && !!U;
	}

	boost::shared_ptr<G3SkyMapWithWeights> Clone(bool copy_data) const {
		if (copy_data)
			return boost::make_shared<G3SkyMapWithWeights>(*this);
		else
			return boost::make_shared<G3SkyMapWithWeights>(this->T,
			    this->IsWeighted(), this->IsPolarized());
	}

	void RemoveWeights();
	void ApplyWeights(G3SkyMapWeightsPtr weights);

	StokesVector get_interp_value(double alpha, double delta) const;

	boost::shared_ptr<G3SkyMapWithWeights> rebin(size_t scale) const;

private:
	template <class A> void serialize(A &ar, const unsigned v);
	friend class cereal::access;
};

G3_POINTERS(G3SkyMapWithWeights);

G3_SERIALIZABLE(G3SkyMap, 2);
G3_SERIALIZABLE(G3SkyMapWeights, 2);
G3_SERIALIZABLE(G3SkyMapWithWeights, 1);

#endif

