#include <pybindings.h>
#include <serialization.h>
#include <container_pybindings.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sstream>

#include <dfmux/DfMuxBuilder.h>

template <class A> void DfMuxBoardSamples::serialize(A &ar, const unsigned v)
{
	using namespace cereal;
	
	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("samples", base_class<std::map<int, DfMuxSamplePtr> >(this));
	ar & make_nvp("nmodules", nmodules);

	if (v > 1) {
		ar & make_nvp("nblocks", nblocks);
		ar & make_nvp("nchannels", nchannels);
	} else {
		nblocks = 1;
		nchannels = 128;
	}
}

G3_SERIALIZABLE_CODE(DfMuxBoardSamples);

std::string DfMuxMetaSample::Summary() const {
	std::ostringstream summ;
	int mods = 0;
	for (auto i = begin(); i != end(); i++)
		mods += i->second.nmodules;
	summ << size() << " boards, with " << mods << " modules";
	return summ.str();
}

template <class A> void DfMuxMetaSample::serialize(A &ar, const unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("boards", base_class<std::map<int, DfMuxBoardSamples> >(this));
}

G3_SERIALIZABLE_CODE(DfMuxMetaSample);

DfMuxBuilder::DfMuxBuilder(int nboards, int64_t collation_tolerance) : 
    G3EventBuilder(MAX_DATASOURCE_QUEUE_SIZE*3), nboards_(nboards),
    last_out_pkt_(0), tolerance_(collation_tolerance), out_pkt_cnt_(0)
{
}

DfMuxBuilder::DfMuxBuilder(std::vector<int> boards,
  int64_t collation_tolerance) : 
    G3EventBuilder(MAX_DATASOURCE_QUEUE_SIZE*3), nboards_(boards.size()),
    last_out_pkt_(0), board_list_(boards), tolerance_(collation_tolerance),
    out_pkt_cnt_(0)
{
}

DfMuxBuilder::~DfMuxBuilder()
{
}

void DfMuxBuilder::ProcessNewData()
{
	DfMuxSamplePacketConstPtr pkt;
	{
		std::lock_guard<std::mutex> lock(queue_lock_);
		pkt = std::dynamic_pointer_cast<const DfMuxSamplePacket>(
		    queue_.front().second);
		queue_.pop_front();
	}
	if (!pkt) {
		log_error("Non-DfMux async data received. Throwing it away.");
		return;
	}

	// Mask data if asked
	if (board_list_.size() > 0 && std::find(board_list_.begin(),
	    board_list_.end(), pkt->board) == board_list_.end())
		return;

	std::deque<struct oqueue_entry>::iterator sample;
	G3TimeStamp timecode;

	// Locate any existing sample at the right time
	timecode = pkt->sample->Timestamp.time;
	for (sample = oqueue_.begin(); sample != oqueue_.end(); sample++) {
		if (llabs(timecode - sample->time) < tolerance_)
			break;
	}

	// Place in queue if at a new time
	if (sample == oqueue_.end()) {
		struct oqueue_entry metasamp;

		// Check if this is some ghost from the past. Log these.
		if (timecode < last_out_pkt_) {
			log_warn("Bogon packet from board %d at past time %s "
			    "(last outbound frame was at %s)", pkt->board,
			    pkt->sample->Timestamp.Description().c_str(),
			    G3Time(last_out_pkt_).Description().c_str());
			return;
		}

		// Figure out where this needs to go in the queue
		for (sample = oqueue_.end(); sample != oqueue_.begin();
		    sample--) {
			if ((sample == oqueue_.begin() ||
			     timecode > (sample-1)->time) &&
			    (sample == oqueue_.end() ||
			     timecode < sample->time))
				break;
		}

		metasamp.frame = std::make_shared<G3Frame>(G3Frame::Timepoint);
		metasamp.sample = std::make_shared<DfMuxMetaSample>();
		metasamp.time = timecode;
		metasamp.frame->Put("EventHeader",
		    std::make_shared<G3Time>(timecode));
		sample = oqueue_.insert(sample, metasamp);
		CollectPolledData(metasamp.frame);
	}

	// Add to meta sample
	int idx = pkt->module * pkt->nblocks + pkt->block;
	if ((*sample->sample)[pkt->board].find(idx) !=
	    (*sample->sample)[pkt->board].end()) {
		log_error("Duplicate packet from board %d module %d/%d at time %s",
		    pkt->board, pkt->module, pkt->block,
		    pkt->sample->Timestamp.Description().c_str());
		return;
	}
	(*sample->sample)[pkt->board][idx] = pkt->sample;
	(*sample->sample)[pkt->board].nmodules = pkt->nmodules;
	(*sample->sample)[pkt->board].nblocks = pkt->nblocks;
	(*sample->sample)[pkt->board].nchannels = pkt->nchannels;
	
	while (oqueue_.size() > 0 && oqueue_.front().sample->size() ==
	    nboards_) {
		std::map<int, DfMuxBoardSamples>::const_iterator i;
		for (i = oqueue_.front().sample->begin();
		    i != oqueue_.front().sample->end(); i++) {
			if (!i->second.Complete())
				break;
		}
		
		if (i == oqueue_.front().sample->end()) {
			oqueue_.front().frame->Put("DfMux",
			    oqueue_.front().sample);
			last_out_pkt_ = oqueue_.front().frame->
			    Get<G3Time>("EventHeader")->time;
			FrameOut(oqueue_.front().frame);
			oqueue_.pop_front();
			out_pkt_cnt_++;
			continue;
		}

		break;
	}

	if (oqueue_.size() >= MAX_DATASOURCE_QUEUE_SIZE) {
		if (out_pkt_cnt_ > 1) {
			// Silently drop the first incomplete sample, which
			// is likely just this code being misaligned with the
			// first packet set. Warn about and save all subsequent
			// incomplete samples.

			// Prepare message about which boards are missing
			// if we have a board list.
			std::ostringstream msg_boards_missing;
			size_t boards_received = nboards_;
			if (board_list_.size() > 0)
				msg_boards_missing << ". Missing boards";

			for (auto b = board_list_.begin();
			    b != board_list_.end(); b++) {
				if (oqueue_.front().sample->find(*b) == 
				    oqueue_.front().sample->end()) {
					msg_boards_missing << " " << *b;
					boards_received -= 1;
					// Insert blank sample from this board
					// so later code knows what could be
					// there.
					(void)(*oqueue_.front().sample)[*b];
				}
			} 
			log_warn("Abandoning missing packets (%s, data "
			    "from %zd/%zd boards received)%s",
			    oqueue_.front().frame->Get<G3Time>(
			      "EventHeader")->Description().c_str(),
			    boards_received, nboards_,
			    msg_boards_missing.str().c_str());

			oqueue_.front().frame->Put("DfMux",
			    oqueue_.front().sample);
			last_out_pkt_ = oqueue_.front().frame->
			    Get<G3Time>("EventHeader")->time;
			FrameOut(oqueue_.front().frame);
		}
		oqueue_.pop_front();
		out_pkt_cnt_++;
	}
}

PYBINDINGS("dfmux", scope)
{
	register_g3map<DfMuxBoardSamples>(scope, "DfMuxBoardSamples",
	  "Container structure for samples from modules on one board, mapping "
	  "0-indexed module and block IDs to a dfmux.DfMuxSample.")
	    .def_readwrite("nmodules", &DfMuxBoardSamples::nmodules,
	      "Number of modules expected to report from this board")
	    .def_readwrite("nblocks", &DfMuxBoardSamples::nblocks,
	      "Number of sub-module blocks expected to report from this board")
	    .def_readwrite("nchannels", &DfMuxBoardSamples::nchannels,
	      "Number of channels per block expected to report from this board")
	    .def("Complete", &DfMuxBoardSamples::Complete,
	      "True if this structure contains data from all expected modules and blocks")
	;

	register_g3map<DfMuxMetaSample>(scope, "DfMuxMetaSample",
	  "Container structure for coincident samples from all boards. "
	  "Individual board data, stored in dfmux.DfMuxBoardSamples classes, "
	  "is contained indexed by board serial number.")
	;

	register_g3module<DfMuxBuilder, G3EventBuilder>(scope, "DfMuxBuilder",
	  "Processing module for data from DfMux boards. Reads data from "
	  "boards data acquisition boards, requiring that data from all "
	  "be timestamped to within collation_tolerance (default 10 "
	  "microseconds) to be considered part of a single sample. If boards "
	  "is an integer, listens for that number. If a list of integers, "
	  "DfMuxBuilder will filter for only boards with serial numbers "
	  "in the list.")
	    .def(py::init<int, int64_t>(), py::arg("boards"),
	        py::arg("collation_tolerance")=10*G3Units::ms)
	    .def(py::init<std::vector<int>, int64_t>(), py::arg("boards"),
	        py::arg("collation_tolerance")=10*G3Units::ms)
	;
}

