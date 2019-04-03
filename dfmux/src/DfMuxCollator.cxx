#include <pybindings.h>
#include <G3Module.h>
#include <G3Data.h>
#include <G3Timestream.h>

#include <dfmux/DfMuxBuilder.h>
#include <dfmux/HardwareMap.h>

#include <sys/types.h>
#include <math.h>

class DfMuxCollator : public G3Module {
public:
	DfMuxCollator(bool flac_compress = true, bool drop_timepoints = true,
	    bool record_sample_times = true);

	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);

private:
	SET_LOGGER("DfMuxCollator");

	G3FramePtr scan_;
	std::deque<G3FramePtr> stash_;
	DfMuxWiringMapConstPtr hwm_;
	bool flac_compress_;
	bool drop_timepoints_;
	bool record_sample_times_;
};

DfMuxCollator::DfMuxCollator(bool flac_compress, bool drop_timepoints,
  bool record_sample_times) :
    G3Module(), flac_compress_(flac_compress),
    drop_timepoints_(drop_timepoints), record_sample_times_(record_sample_times)
{
}

struct bolots_cache_item {
	G3TimestreamPtr i;
	G3TimestreamPtr q;
	const DfMuxChannelMapping *mapping;
};

void DfMuxCollator::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
{
	if (frame->type == G3Frame::Wiring) {
		hwm_ = frame->Get<DfMuxWiringMap>("WiringMap");
		out.push_back(frame);
		return;
	} else if (frame->type == G3Frame::Timepoint) {
		// Store timepoints iff the first scan frame has come by
		if (scan_)
			stash_.push_back(frame);
		else if (!drop_timepoints_)
			out.push_back(frame);
		return;
	} else if (frame->type != G3Frame::Scan &&
	    frame->type != G3Frame::EndProcessing) {
		out.push_back(frame);
		return;
	}

	// If we get an EndProcessing frame immediately with no in-progress
	// scan, just bail
	if (frame->type == G3Frame::EndProcessing && !scan_) {
		out.push_back(frame);
		return;
	}

	// If we got this far, the frame is either a new scan or an end-of-file
	// signal. Begin collating in either case.
	G3TimestreamMapPtr itsmap(new G3TimestreamMap),
	    qtsmap(new G3TimestreamMap);
	G3VectorTimePtr sample_times(new G3VectorTime);
	std::map<std::string, G3TimestreamPtr> extra_data;
	G3Time start, stop;
	int sample = 0;

	if (!hwm_)
		log_fatal("No hardware wiring map so far");

	// Init bolo timestreams and make a small iterable cache to reduce
	// number of find() calls below.
	G3Timestream ts_base(stash_.size(), NAN);
	ts_base.units = G3Timestream::Counts;
	if (flac_compress_)
		ts_base.SetFLACCompression(5);
	std::vector<bolots_cache_item> hwmts_cache;
	hwmts_cache.reserve(hwm_->size());
	for (auto chan = hwm_->begin(); chan != hwm_->end(); chan++) {
		bolots_cache_item it;
		it.i = (*itsmap)[chan->first] = G3TimestreamPtr(
		    new G3Timestream(ts_base));
		it.q = (*qtsmap)[chan->first] = G3TimestreamPtr(
		    new G3Timestream(ts_base));
		it.mapping = &chan->second;
		hwmts_cache.push_back(it);
	}

	// For the first scan, don't process the nonexistant previous scan
	if (!scan_ || stash_.size() == 0)
		goto out;

	// Init aux timestreams by enumerating all G3Doubles in first point
	{
		std::vector<std::string> keys = (*stash_.begin())->Keys();
		for (auto key = keys.begin(); key != keys.end(); key++) {
			// Make a timestream for each. Note that there are no
			// guarantees about bit width, so don't turn on FLAC
			// here.
			if ((*stash_.begin())->Has<G3Double>(*key))
				extra_data[*key] =
				    boost::make_shared<G3Timestream>(
				    stash_.size(), NAN);
		}
	}

	// Now collect all the data
	start = *((*stash_.begin())->Get<G3Time>("EventHeader"));
	for (auto t = stash_.begin(); t != stash_.end(); t++, sample++) {
		DfMuxMetaSampleConstPtr metasamp =
		    (*t)->Get<DfMuxMetaSample>("DfMux");

		// Collect timestamps
		stop = *((*t)->Get<G3Time>("EventHeader"));
		sample_times->push_back(stop);

		// Skip points that have been marked as garbage (e.g. by the
		// scanifier).  Implicitly fills these points with NaNs.
		if ((*t)->Has("GarbageData")) {
			if ((*t)->Get<G3Bool>("GarbageData")) {
				t->reset();
				continue;
			}
		}

		// Collect all the bolo data
		for (auto chan = hwmts_cache.begin(); chan != hwmts_cache.end();
		    chan++) {
			auto board = metasamp->find(
			    chan->mapping->board_serial);
			if (board == metasamp->end())
				continue;

			auto module = board->second.find(chan->mapping->module);
			if (module == board->second.end())
				continue;

			if (module->second->size()/2 <
			    size_t(chan->mapping->channel)) {
				log_fatal("Board %d, module %d only has %zd "
				    "channels, but trying to read %d",
				    chan->mapping->board_serial,
				    chan->mapping->module,
				    module->second->size()/2,
				    chan->mapping->channel);
			}

			(*chan->i)[sample] =
			    (*module->second)[chan->mapping->channel*2];
			(*chan->q)[sample] =
			    (*module->second)[chan->mapping->channel*2 + 1];
		}

		// Next run through aux data. Missing points get filled in
		// implicitly with NaN, as above.
		for (auto aux = extra_data.begin(); aux != extra_data.end();
		    aux++) {
			if ((*t)->Has<G3Double>(aux->first))
				(*aux->second)[sample] =
				    (*t)->Get<G3Double>(aux->first)->value;
		}

		// If dropping the timepoint frames anyway, free memory now
		t->reset();
	}

	// Set begin and end times
	for (auto chan = hwmts_cache.begin(); chan != hwmts_cache.end();
	    chan++) {
		chan->i->start = chan->q->start = start;
		chan->i->stop = chan->q->stop = stop;
	}
	for (auto aux = extra_data.begin(); aux != extra_data.end(); aux++) {
		aux->second->start = start;
		aux->second->stop = stop;
	}

	if (drop_timepoints_)
		stash_.clear();
	else
		stash_.swap(out);

	// Load bolo and aux data
	scan_->Put("RawTimestreams_I", itsmap);
	scan_->Put("RawTimestreams_Q", qtsmap);

	for (auto aux = extra_data.begin(); aux != extra_data.end(); aux++)
		scan_->Put(aux->first, aux->second);

	// Save timestamps if requested
	if (record_sample_times_)
		scan_->Put("DetectorSampleTimes", sample_times);

	out.push_back(scan_);

out:
	if (frame->type == G3Frame::Scan)
		scan_ = frame;

	if (frame->type == G3Frame::EndProcessing)
		out.push_back(frame);
}

EXPORT_G3MODULE("dfmux", DfMuxCollator, (init<optional<bool, bool, bool> >(
    args("flac_compress", "drop_timepoints", "record_sampletimes"))),
    "Collects DfMux timepoints into scan frames using a provided wiring map. "
    "Scan frames are created when an empty Scan frame appears in the data "
    "stream. This frame will contain all subsequent timepoints until either "
    "the next Scan frame is detected or the stream ends. In addition to dfmux "
    "timestreams, any scalar floating numbers that recur in every input "
    "Timepoint frame will be combined into a G3Timestream of the same name "
    "stored in the output scan frame."
);

