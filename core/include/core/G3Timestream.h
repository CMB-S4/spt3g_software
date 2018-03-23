#ifndef _G3_TIMESTREAM_H
#define _G3_TIMESTREAM_H

#include <G3Frame.h>
#include <G3TimeStamp.h>
#include <G3Vector.h>
#include <vector>
#include <map>

class G3Timestream : public G3Vector<double> {
public:
	enum TimestreamUnits {
		None = 0,
		Counts = 1,
		Current = 2,
		Power = 3,
		Resistance = 5,
		Tcmb = 4,
	};

	G3Timestream() : G3Vector<double>(), units(None), use_flac_(0) {}
	G3Timestream(std::vector<double>::size_type s) : G3Vector<double>(s),
	    units(None), use_flac_(0) {}
	G3Timestream(std::vector<double>::size_type s,
	    const double &val) : G3Vector<double>(s, val), units(None),
	    use_flac_(0) {}
	G3Timestream(const G3Timestream &r) : G3Vector<double>(r),
	    units(r.units), start(r.start), stop(r.stop),
	    use_flac_(r.use_flac_) {}
	template <typename Iterator> G3Timestream(Iterator l, Iterator r) :
	    G3Vector<double>(l, r), units(None), use_flac_(0) {}

	// FLAC compression levels range from 0-9. 0 means do not use FLAC.
	void SetFLACCompression(int compression_level);

	TimestreamUnits units;
	G3Time start, stop;

	double GetSampleRate() const;
	
	template <class A> void load(A &ar, unsigned v);
	template <class A> void save(A &ar, unsigned v) const;

	std::string Description() const;
	std::string Summary() const { return Description(); };

	// Arithmetic operations
	G3Timestream operator +(const G3Timestream &r) const;
	G3Timestream &operator +=(const G3Timestream &r);
	G3Timestream &operator -=(const G3Timestream &r);
	G3Timestream operator -(const G3Timestream &r) const;
	G3Timestream operator *(const G3Timestream &r) const;
	G3Timestream operator /(const G3Timestream &r) const;
	G3Timestream operator +(double r) const;
	G3Timestream operator -(double r) const;
	G3Timestream operator *(double r) const;
	G3Timestream &operator *=(double r);
	G3Timestream operator /(double r) const;
	G3Timestream &operator /=(double r);
private:
	SET_LOGGER("G3Timestream");

	uint8_t use_flac_;
};

G3_POINTERS(G3Timestream);

class G3TimestreamMap : public G3FrameObject,
    public std::map<std::string, G3TimestreamPtr> {
public:
	// Return true if all timestreams start and end at the same time and
	// contain the same number of samples.
	bool CheckAlignment() const;

	// Return statistics for contained timestreams. Results undefined if
	// CheckAlignment is false.
	G3Time GetStartTime() const;
	G3Time GetStopTime() const;
	double GetSampleRate() const;
	size_t NSamples() const;

	template <class A> void serialize(A &ar, unsigned v);
	std::string Description() const;
};

G3_POINTERS(G3TimestreamMap);

namespace cereal {
	template <class A> struct specialize<A, G3Timestream, cereal::specialization::member_load_save> {};
	template <class A> struct specialize<A, G3TimestreamMap, cereal::specialization::member_serialize> {};
}

G3_SERIALIZABLE(G3Timestream, 2);
G3_SERIALIZABLE(G3TimestreamMap, 3);

#endif

