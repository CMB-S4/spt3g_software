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

	G3SkyMap(MapCoordReference coords, bool weighted_ = true,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    MapPolType pol_type = None) :
	    coord_ref(coords), units(u), pol_type(pol_type),
	    weighted(weighted_), overflow(0) {}
	virtual ~G3SkyMap() {};

	// Reimplement the following in subclasses
	virtual boost::shared_ptr<G3SkyMap> Clone(bool copy_data = true) const = 0;

	MapCoordReference coord_ref;
	G3Timestream::TimestreamUnits units;
	MapPolType pol_type;
	bool weighted;
	double overflow;

	void SetOverflow(double val) __attribute__((deprecated)) {
		overflow = val;
	}
	double GetOverflow(void) const __attribute__((deprecated)) {
		return overflow;
	}

	// Return a (modifiable) pixel value
	virtual double &operator [] (size_t i) = 0;
	double operator [] (size_t i) const { return this->at(i); };
	virtual double at(size_t i) const = 0;

	virtual size_t size(void) const;  // total number of pixels
	virtual size_t NpixAllocated(void) const = 0;  // stored in RAM
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
	std::vector<int> AnglesToPixels(const std::vector<double> & alphas,
	    const std::vector<double> & deltas) const;
	void PixelsToAngles(const std::vector<int> & pixels,
	    std::vector<double> & alphas, std::vector<double> & deltas) const;

	virtual std::vector<double> PixelToAngle(size_t pixel) const = 0;
	virtual size_t AngleToPixel(double alpha, double delta) const = 0;

	// Rebinning and interpolation
	virtual void GetRebinAngles(long pixel, size_t scale,
	    std::vector<double> & alphas, std::vector<double> & deltas) const = 0;
	virtual void GetInterpPixelsWeights(double alpha, double delta,
	    std::vector<long> & pixels, std::vector<double> & weights) const = 0;
	double GetInterpPrecalc(const std::vector<long> &pixels,
	    const std::vector<double> &weights) const;
	double GetInterpValue(double alpha, double delta) const;
	std::vector<double> GetInterpValues(const std::vector<double> &alphas,
	    const std::vector<double> &deltas) const;

	virtual boost::shared_ptr<G3SkyMap> Rebin(size_t scale, bool norm = true) const = 0;

	virtual bool IsDense() const {
		throw std::runtime_error("Checking array density not implemented");
	}
	virtual void ConvertToDense() {
		throw std::runtime_error("Conversion to dense array not implemented");
	}
	virtual void Compress(bool zero_nans = false) {
		throw std::runtime_error("Compression not implemented");
	}

protected:
	virtual void InitFromV1Data(std::vector<size_t>,
	    const std::vector<double> &) {
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
	StokesVector(double &t_) : t(t_), q(backing[1]), u(backing[2]) {}
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

	StokesVector &operator /=(const double r) {
		t /= r; q /= r; u /= r;
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
	MuellerMatrix(double &tt_) : tt(tt_), tq(backing[1]), tu(backing[2]),
	    qq(backing[3]), qu(backing[4]), uu(backing[5]) {}
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

	MuellerMatrix &operator /=(const double r) {
		tt /= r; tq /= r; tu /= r;
		qq /= r; qu /= r; uu /= r;
		return *this;
	}

	StokesVector operator *(const StokesVector &r) const {
		StokesVector s;
		s.t = tt * r.t + tq * r.q + tu * r.u;
		s.q = tq * r.t + qq * r.q + qu * r.u;
		s.u = tu * r.t + qu * r.q + uu * r.u;
		return s;
	}

	double Det() const {
		return (tt * (qq * uu - qu * qu) -
			tq * (tq * uu - qu * tu) +
			tu * (tq * qu - qq * tu));
	}

	MuellerMatrix Inv() const;
	double Cond() const;

private:
	double backing[6];
};

class G3SkyMapWeights : public G3FrameObject {
public:
	G3SkyMapWeights() {}

	// Instantiate weight maps based on the metadata of a reference map
	G3SkyMapWeights(G3SkyMapConstPtr ref_map, bool polarized = true);

	G3SkyMapWeights(const G3SkyMapWeights &r);
	G3SkyMapPtr TT, TQ, TU, QQ, QU, UU;

	bool IsPolarized() const {
		return !!TQ && !!TU && !!QQ && !!QU && !!UU;
	}

	bool IsCongruent() const {
		if (!TT || !IsPolarized())
			return true;
		return (TT->IsCompatible(*TQ) &&
			TT->IsCompatible(*TU) &&
			TT->IsCompatible(*QQ) &&
			TT->IsCompatible(*QU) &&
			TT->IsCompatible(*UU));
	}

	MuellerMatrix operator [] (int i) {
		return (!IsPolarized()) ? MuellerMatrix((*TT)[i]) :
		    MuellerMatrix((*TT)[i], (*TQ)[i], (*TU)[i], (*QQ)[i],
		        (*QU)[i], (*UU)[i]);
	}

	const MuellerMatrix at (int i) const {
		MuellerMatrix m;
		m.tt = TT->at(i);
		if (IsPolarized()) {
			m.tq = TQ->at(i);
			m.tu = TU->at(i);
			m.qq = QQ->at(i);
			m.qu = QU->at(i);
			m.uu = UU->at(i);
		}
		return m;
	}

	size_t size(void) const { return TT->size(); };  // total number of pixels
	std::vector<size_t> shape(void) const { return TT->shape(); };  // map shape

	// Coadd
	G3SkyMapWeights &operator+=(const G3SkyMapWeights &rhs);

	// Scale
	G3SkyMapWeights &operator*=(double rhs);
	G3SkyMapWeights &operator/=(double rhs);

	// Mask
	G3SkyMapWeights &operator*=(const G3SkyMap &rhs);

	G3SkyMapPtr Det() const;
	G3SkyMapPtr Cond() const;

	boost::shared_ptr<G3SkyMapWeights> Rebin(size_t scale) const;

	void Compress(bool zero_nans = false);

	boost::shared_ptr<G3SkyMapWeights> Clone(bool copy_data) const {
		if (copy_data)
			return boost::make_shared<G3SkyMapWeights>(*this);
		else
			return boost::make_shared<G3SkyMapWeights>(this->TT, this->IsPolarized());
	}
private:
	template <class A> void serialize(A &ar, const unsigned v);
	friend class cereal::access;
};

G3_POINTERS(G3SkyMapWeights);

class G3SkyMapWithWeights : public G3FrameObject {
public:
	G3SkyMapWithWeights() {}

	G3SkyMapWithWeights(G3SkyMapConstPtr ref_map, bool weighted = true,
	    bool polarized = true, std::string map_id = "");

	G3SkyMapWithWeights(const G3SkyMapWithWeights &r);

	G3SkyMapPtr T, Q, U;
	G3SkyMapWeightsPtr weights;

	std::string map_id;

	StokesVector operator [] (int i) {
		return IsPolarized() ? StokesVector((*T)[i], (*Q)[i], (*U)[i]) :
		    StokesVector((*T)[i]);
	}

	const StokesVector at (int i) const {
		StokesVector v;
		v.t = T->at(i);
		v.q = !Q ? 0 : Q->at(i);
		v.u = !U ? 0 : U->at(i);
		return v;
	}

	size_t size(void) const { return T->size(); };  // total number of pixels
	std::vector<size_t> shape(void) const { return T->shape(); };  // map shape

	// Coadd
	G3SkyMapWithWeights &operator+=(const G3SkyMapWithWeights &rhs);

	// Scale
	G3SkyMapWithWeights &operator*=(double rhs);
	G3SkyMapWithWeights &operator/=(double rhs);

	// Mask
	G3SkyMapWithWeights &operator*=(const G3SkyMap &rhs);

	bool IsWeighted() const {
		return !!weights && !!weights->TT;
	}

	bool IsPolarized() const {
		return !!Q && !!U;
	}

	bool IsCongruent() const {
		if (!T)
			return true;
		if (IsWeighted()) {
			if  (!(weights->TT))
				return  false;
			if (!(T->IsCompatible(*(weights->TT))))
				return false;
			if (IsPolarized() != weights->IsPolarized())
				return false;
			if (!(weights->IsCongruent()))
				return false;
		}
		if (IsPolarized())
			return T->IsCompatible(*Q) && T->IsCompatible(*U);
		return true;
	}

	boost::shared_ptr<G3SkyMapWithWeights> Clone(bool copy_data) const {
		if (copy_data)
			return boost::make_shared<G3SkyMapWithWeights>(*this);
		else
			return boost::make_shared<G3SkyMapWithWeights>(this->T,
			    this->IsWeighted(), this->IsPolarized());
	}

	G3SkyMapWeightsPtr RemoveWeights();
	void ApplyWeights(G3SkyMapWeightsPtr weights);

	boost::shared_ptr<G3SkyMapWithWeights> Rebin(size_t scale) const;

	void Compress(bool zero_nans = false);

private:
	template <class A> void serialize(A &ar, const unsigned v);
	friend class cereal::access;
};

G3_POINTERS(G3SkyMapWithWeights);

G3_SERIALIZABLE(G3SkyMap, 2);
G3_SERIALIZABLE(G3SkyMapWeights, 3);
G3_SERIALIZABLE(G3SkyMapWithWeights, 1);

#endif

