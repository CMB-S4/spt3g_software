#ifndef _G3_SKYMAP_H
#define _G3_SKYMAP_H

#include <G3Frame.h>
#include <G3Timestream.h>
#include <G3Quat.h>
#include <vector>
#include <map>

#include <maps/G3SkyMapMask.h>

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

#define NOT_IMPLEMENTED { log_fatal("Not implemented"); }

class G3SkyMap {
public:
	// Following numerical values are important
	enum MapPolType {
		T = 0,
		Q = 1,
		U = 2,
		E = 3,
		B = 4,
		None = 7,
		TT = 8,
		TQ = 9,
		TU = 10,
		QQ = 11,
		QU = 12,
		UU = 13,
		TE = 14,
		TB = 15,
		EE = 16,
		EB = 17,
		BB = 18,
	};

	enum MapPolConv {
		IAU = 0,
		COSMO = 1,
		ConvNone = 2
	};

	G3SkyMap(MapCoordReference coords, bool weighted_ = true,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    MapPolType pol_type = None, MapPolConv pol_conv = ConvNone);
	virtual ~G3SkyMap() {};

	// Reimplement the following in subclasses
	virtual std::shared_ptr<G3SkyMap> Clone(bool copy_data = true) const NOT_IMPLEMENTED;

	MapCoordReference coord_ref;
	G3Timestream::TimestreamUnits units;
	MapPolType pol_type;
	MapPolConv pol_conv;
	bool weighted;
	double overflow;

	void SetOverflow(double val) __attribute__((deprecated)) {
		overflow = val;
	}
	double GetOverflow(void) const __attribute__((deprecated)) {
		return overflow;
	}

	// Return a (modifiable) pixel value
	virtual double &operator [] (size_t i) NOT_IMPLEMENTED;
	double operator [] (size_t i) const { return this->at(i); };
	virtual double at(size_t i) const NOT_IMPLEMENTED;

	virtual size_t size(void) const;  // total number of pixels
	virtual size_t NpixAllocated(void) const NOT_IMPLEMENTED;  // stored in RAM
	virtual size_t NpixNonZero(void) const NOT_IMPLEMENTED;  // nonzero
	virtual std::vector<size_t> shape(void) const NOT_IMPLEMENTED;  // map shape

	bool IsPolarized() const { return pol_conv != ConvNone; }

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
	virtual G3SkyMap &operator*=(const G3SkyMapMask &rhs);
	virtual G3SkyMap &operator*=(double rhs);

	// /
	virtual G3SkyMap &operator/=(const G3SkyMap &rhs);
	virtual G3SkyMap &operator/=(double rhs);

	// Comparison operations

	// <
	virtual G3SkyMapMask operator<(const G3SkyMap &rhs) const;
	virtual G3SkyMapMask operator<=(const G3SkyMap &rhs) const;
	virtual G3SkyMapMask operator<(double rhs) const;
	virtual G3SkyMapMask operator<=(double rhs) const;

	// ==
	virtual G3SkyMapMask operator==(const G3SkyMap &rhs) const;
	virtual G3SkyMapMask operator!=(const G3SkyMap &rhs) const;
	virtual G3SkyMapMask operator==(double rhs) const;
	virtual G3SkyMapMask operator!=(double rhs) const;

	// <
	virtual G3SkyMapMask operator>(const G3SkyMap &rhs) const;
	virtual G3SkyMapMask operator>=(const G3SkyMap &rhs) const;
	virtual G3SkyMapMask operator>(double rhs) const;
	virtual G3SkyMapMask operator>=(double rhs) const;

	virtual bool all(G3SkyMapMaskConstPtr where=NULL) const;
	virtual bool any(G3SkyMapMaskConstPtr where=NULL) const;
	virtual double sum(G3SkyMapMaskConstPtr where=NULL) const;
	virtual double mean(G3SkyMapMaskConstPtr where=NULL) const;
	virtual double median(G3SkyMapMaskConstPtr where=NULL) const;
	virtual double var(size_t ddof, G3SkyMapMaskConstPtr where=NULL) const;
	virtual double std(size_t ddof, G3SkyMapMaskConstPtr where=NULL) const;
	virtual double min(G3SkyMapMaskConstPtr where=NULL) const;
	virtual double max(G3SkyMapMaskConstPtr where=NULL) const;
	virtual size_t argmin(G3SkyMapMaskConstPtr where=NULL) const;
	virtual size_t argmax(G3SkyMapMaskConstPtr where=NULL) const;

	virtual double nansum(G3SkyMapMaskConstPtr where=NULL) const;
	virtual double nanmean(G3SkyMapMaskConstPtr where=NULL) const;
	virtual double nanmedian(G3SkyMapMaskConstPtr where=NULL) const;
	virtual double nanvar(size_t ddof, G3SkyMapMaskConstPtr where=NULL) const;
	virtual double nanstd(size_t ddof, G3SkyMapMaskConstPtr where=NULL) const;
	virtual double nanmin(G3SkyMapMaskConstPtr where=NULL) const;
	virtual double nanmax(G3SkyMapMaskConstPtr where=NULL) const;
	virtual size_t nanargmin(G3SkyMapMaskConstPtr where=NULL) const;
	virtual size_t nanargmax(G3SkyMapMaskConstPtr where=NULL) const;

	virtual G3SkyMapMask isinf(G3SkyMapMaskConstPtr where=NULL) const;
	virtual G3SkyMapMask isnan(G3SkyMapMaskConstPtr where=NULL) const;
	virtual G3SkyMapMask isfinite(G3SkyMapMaskConstPtr where=NULL) const;

	virtual void ApplyMask(const G3SkyMapMask &mask, bool inverse=false);

	// Pointing information
	std::vector<uint64_t> AnglesToPixels(const std::vector<double> & alphas,
	    const std::vector<double> & deltas) const;
	void PixelsToAngles(const std::vector<uint64_t> & pixels,
	    std::vector<double> & alphas, std::vector<double> & deltas) const;
	std::vector<uint64_t> QuatsToPixels(const G3VectorQuat &quats) const;
	G3VectorQuat PixelsToQuats(const std::vector<uint64_t> &pixels) const;

	virtual std::vector<double> PixelToAngle(size_t pixel) const;
	virtual size_t AngleToPixel(double alpha, double delta) const;
	virtual Quat PixelToQuat(size_t pixel) const NOT_IMPLEMENTED;
	virtual size_t QuatToPixel(const Quat &q) const NOT_IMPLEMENTED;

	// Rebinning and interpolation
	virtual void GetRebinAngles(size_t pixel, size_t scale,
	    std::vector<double> & alphas, std::vector<double> & deltas) const;
	virtual G3VectorQuat GetRebinQuats(size_t pixel, size_t scale) const NOT_IMPLEMENTED;
	virtual void GetInterpPixelsWeights(double alpha, double delta,
	    std::vector<uint64_t> & pixels, std::vector<double> & weights) const;
	virtual void GetInterpPixelsWeights(const Quat &q, std::vector<uint64_t> & pixels,
	    std::vector<double> & weights) const NOT_IMPLEMENTED;
	double GetInterpPrecalc(const std::vector<uint64_t> &pixels,
	    const std::vector<double> &weights) const;
	double GetInterpValue(double alpha, double delta) const;
	double GetInterpValue(const Quat &q) const;
	std::vector<double> GetInterpValues(const std::vector<double> &alphas,
	    const std::vector<double> &deltas) const;
	std::vector<double> GetInterpValues(const G3VectorQuat & quats) const;

	virtual std::shared_ptr<G3SkyMap> Rebin(size_t scale, bool norm = true) const NOT_IMPLEMENTED;

	/* Analogue to healpy.query_disc, returns list of pixels within a disc */
	std::vector<uint64_t> QueryDisc(double alpha, double delta, double radius) const;
	virtual std::vector<uint64_t> QueryDisc(const Quat &q, double radius) const NOT_IMPLEMENTED;
	std::vector<uint64_t> QueryAlphaEllipse(double alpha, double delta, double a, double b) const;
	std::vector<uint64_t> QueryAlphaEllipse(const Quat &q, double a, double b) const;

	virtual bool IsDense() const {
		throw std::runtime_error("Checking array density not implemented");
	}
	virtual void ConvertToDense() {
		throw std::runtime_error("Conversion to dense array not implemented");
	}
	virtual void Compact(bool zero_nans = false) {
		throw std::runtime_error("Compactification not implemented");
	}

	virtual G3SkyMapMaskPtr MakeMask(bool zero_nans = false,
	    bool zero_infs = false) const;

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

#undef NOT_IMPLEMENTED

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

	// Constructor for polarization coupling
	StokesVector(double pol_ang, double pol_eff);

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

	void rotate_pol(double pol_ang);

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

	// Constructor from Stokes vector
	MuellerMatrix(const StokesVector &r) : tt(backing[0]), tq(backing[1]),
	    tu(backing[2]), qq(backing[3]), qu(backing[4]), uu(backing[5]) {
		from_vector(r);
	}

	double &tt, &tq, &tu, &qq, &qu, &uu;

	void from_vector(const StokesVector &r) {
		tt = r.t * r.t; tq = r.t * r.q; tu = r.t * r.u;
		qq = r.q * r.q; qu = r.q * r.u; uu = r.u * r.u;
	}

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
	MuellerMatrix operator *(const double r) {
		MuellerMatrix ret; // No copy to initialize to backing
		ret.tt = tt*r; ret.tq = tq*r; ret.tu = tu*r;
		ret.qq = qq*r; ret.qu = qu*r; ret.uu = uu*r;
		return ret;
	}

	MuellerMatrix &operator /=(const double r) {
		tt /= r; tq /= r; tu /= r;
		qq /= r; qu /= r; uu /= r;
		return *this;
	}
	MuellerMatrix operator /(const double r) {
		MuellerMatrix ret; // No copy to initialize to backing
		ret.tt = tt/r; ret.tq = tq/r; ret.tu = tu/r;
		ret.qq = qq/r; ret.qu = qu/r; ret.uu = uu/r;
		return ret;
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

	MuellerMatrix inv() const;
	double cond() const;

	void rotate_pol(double pol_ang);

private:
	double backing[6];
};

class G3SkyMapWeights : public G3FrameObject {
public:
	G3SkyMapWeights() {}

	// Instantiate weight maps based on the metadata of a reference map
	G3SkyMapWeights(G3SkyMapConstPtr ref_map);

	G3SkyMapWeights(const G3SkyMapWeights &r, bool copy_data = true);
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

	bool IsCompatible(const G3SkyMap & other) const {
		if (!TT)
			return false;
		return TT->IsCompatible(other);
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
	G3SkyMapWeights &operator-=(const G3SkyMapWeights &rhs);

	// Scale
	G3SkyMapWeights &operator*=(double rhs);
	G3SkyMapWeights &operator/=(double rhs);

	// Mask
	G3SkyMapWeights &operator*=(const G3SkyMap &rhs);
	G3SkyMapWeights &operator*=(const G3SkyMapMask &rhs);

	G3SkyMapPtr Det() const;
	G3SkyMapPtr Cond() const;
	std::shared_ptr<G3SkyMapWeights> Inv() const;

	std::shared_ptr<G3SkyMapWeights> Rebin(size_t scale) const;

	void ApplyMask(const G3SkyMapMask &mask, bool inverse=false);

	void Compact(bool zero_nans = false);

	std::shared_ptr<G3SkyMapWeights> Clone(bool copy_data) const {
		if (copy_data)
			return std::make_shared<G3SkyMapWeights>(*this, true);
		else
			return std::make_shared<G3SkyMapWeights>(*this, false);
	}
private:
	template <class A> void serialize(A &ar, const unsigned v);
	friend class cereal::access;
};

G3_POINTERS(G3SkyMapWeights);

G3_SERIALIZABLE(G3SkyMap, 3);
G3_SERIALIZABLE(G3SkyMapWeights, 3);

#endif

