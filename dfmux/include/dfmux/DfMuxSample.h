#ifndef DFMUXSAMPLE_H
#define DFMUXSAMPLE_H

#include <G3Frame.h>
#include <G3TimeStamp.h>

#include <vector>

/*
 * Class representing a single DfMux sample.  In the context of
 * DfMuxMonitor, a sample represents a single timestamped sample
 * for each channel on each board in the specified system.
 * Samples are (I, Q) pairs for each channel, so that sample 0 is
 * the first channel's I value, 1 is the first channel's Q value,
 * 2 is the second channel's I value, etc.
 */
class DfMuxSample : public G3FrameObject, public std::vector<int32_t>
{
public:
	DfMuxSample() : G3FrameObject(), Timestamp(0) {}
	DfMuxSample(G3TimeStamp timestamp, int nsamples) : G3FrameObject(),
	    std::vector<int32_t>(nsamples), Timestamp(timestamp) {}

	int32_t *Samples() const {return (int32_t *)&(*this)[0];}
	const int NSamples() const {return size();}

	G3Time Timestamp;

	template <class A> void serialize(A &ar, unsigned v);
};

namespace cereal {
	template <class A> struct specialize<A, DfMuxSample, cereal::specialization::member_serialize> {};
}

G3_POINTERS(DfMuxSample);
G3_SERIALIZABLE(DfMuxSample, 1);

#endif

