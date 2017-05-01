#ifndef DFMUX_DATA_BUILDER_H
#define DFMUX_DATA_BUILDER_H

#include <deque>
#include <unordered_map>
#include <G3EventBuilder.h>

#include <dfmux/DfMuxCollector.h>

#define MAX_DATASOURCE_QUEUE_SIZE 1000 // Maximum number of samples to wait before dropping data

/*
 * DfMuxBoardSamples: samples from a board, indexed by module ID (0-7)
 */
struct DfMuxBoardSamples : public G3FrameObject, public std::map<int32_t, DfMuxSamplePtr> {
	size_t nmodules; // Total number of modules expected from this board
	
	bool Complete() const { return (size() == nmodules); };
	template <class A> void serialize(A &ar, unsigned v);
};

G3_POINTERS(DfMuxBoardSamples);
G3_SERIALIZABLE(DfMuxBoardSamples, 1);

/*
 * DfMuxMetaSample: collection of DfMuxBoardSamples, indexed by board serial
 * number. Standard class stored in frames.
 */
struct DfMuxMetaSample : public G3FrameObject, public std::map<int32_t, DfMuxBoardSamples> {
	
	virtual std::string Summary() const;
	template <class A> void serialize(A &ar, unsigned v);
};

G3_POINTERS(DfMuxMetaSample);
G3_SERIALIZABLE(DfMuxMetaSample, 1);

namespace cereal {
	template <class A> struct specialize<A, DfMuxBoardSamples, cereal::specialization::member_serialize> {};
	template <class A> struct specialize<A, DfMuxMetaSample, cereal::specialization::member_serialize> {};
}

/*
 * G3EventBuilder subclass for processing DfMux data. Collects data from
 * DfMuxCollector (or LegacyDfMuxCollector) and organizes it into Timepoint
 * frames containing all the samples from boards within collation_tolerance.
 *
 * Output frames have two keys:
 * - EventHeader: G3Time with the time of the first board sample in the frame.
 * - DfMux: A DfMuxMetaSample with the data from all the boards within the
 *   collation window.
 */
class DfMuxBuilder : public G3EventBuilder {
public:
	DfMuxBuilder(int boards,
	    int64_t collation_tolerance=10*G3Units::microsecond);
	DfMuxBuilder(std::vector<int> boards,
	    int64_t collation_tolerance=10*G3Units::microsecond);

	~DfMuxBuilder();

protected:
	void ProcessNewData();

private:
	struct oqueue_entry {
		G3FramePtr frame;
		DfMuxMetaSamplePtr sample;
		G3TimeStamp time;
	};
	std::deque<struct oqueue_entry> oqueue_;
	std::map<int32_t, std::map<int32_t, int32_t> > sequence_;

	size_t nboards_;
	G3TimeStamp last_out_pkt_;
	std::vector<int> board_list_;
	int64_t tolerance_;
	uint64_t out_pkt_cnt_;

	SET_LOGGER("DfMuxBuilder");
};

G3_POINTERS(DfMuxBuilder);

#endif // DFMUX_DATA_COLLECTOR_H

