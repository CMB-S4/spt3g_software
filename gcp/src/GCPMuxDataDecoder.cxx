#include <pybindings.h>

#include <G3Module.h>
#include <G3TimeStamp.h>
#include <G3Data.h>
#include <G3Map.h>
#include <G3Vector.h>

#include <dfmux/HardwareMap.h>
#include <dfmux/DfMuxSample.h>
#include <dfmux/DfMuxBuilder.h>

#include <boost/tokenizer.hpp>
#include <arpa/inet.h>

// Construct a fake board IP 192.168.1.X. This is used for internal
// bookkeeping and need not be the real IP so long as it is unique.
#define FAKE_BOARD_IP(serial) htonl((192 << 24) | (168 << 16) | \
		    (1 << 8) | (serial))

class GCPMuxDataDecoder : public G3Module {
public:
	GCPMuxDataDecoder() : G3Module(), hwm_emitted_(false) {}

	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);

private:
	G3FramePtr EmitWiringMap(G3FramePtr input);
	DfMuxWiringMapConstPtr cached_wiring_map_;
	
	bool hwm_emitted_;
	size_t max_channel_count_;
	std::vector<DfMuxChannelMapping> channels_;
	std::map<int, int> board_id_map_;

	SET_LOGGER("GCPMuxDataDecoder");
};

EXPORT_G3MODULE("gcp", GCPMuxDataDecoder, init<>(),
    "Extracts contents of receiver registers in SPTpol-style ARC files into a "
    "wiring map and timepoint frames. This is designed to convert SPTpol-style "
    "data in which GCP records bolometer data into the ARC files into a format "
    "equivalent to that for SPT-3G. \n\n"
    "For old (100d SPTpol) data not containing wiring information, insert a "
    "wiring map into the pipeline ahead of this module. The board_serial "
    "should be set to a real (positive) value for all bolometer channels and "
    "-1 for the calibrator sync readout.");

G3FramePtr
GCPMuxDataDecoder::EmitWiringMap(G3FramePtr input)
{
	/*
	 * Transforms receiver.bolometer.dfml_addr and
	 * receiver.bolometer.mbi_addr registers into a WiringMap.
	 * The results are then cached for construction of Timepoint
	 * frames in the main Process() method.
	 *
	 * Also emits a mapping called PhysicalBoloIDs into the
	 * Wiring frame. This corresponds to the physical_name property
	 * of BolometerProperties, which is a calibration quantity, but
	 * GCP does not store the rest of that data. As such, we just
	 * store it and let another module patch things up later.
	 */

	G3MapFrameObjectConstPtr recv, bolos;
	G3VectorFrameObjectConstPtr dfml_addr, mbi_addr, phys_id_reg;

	DfMuxWiringMapPtr wiring(new DfMuxWiringMap);
	G3MapStringPtr phys_ids(new G3MapString);

	if (!input->Has("receiver"))
		log_fatal("Missing receiver registers in ARC file. Is this "
		    "SPT3G data?");

	recv = input->Get<G3MapFrameObject>("receiver");
	bolos = boost::dynamic_pointer_cast<G3MapFrameObject>(
	    recv->at("bolometers"));
	phys_id_reg = boost::dynamic_pointer_cast<G3VectorFrameObject>(
	    bolos->at("id"));
	max_channel_count_ = 0;

	// Get board indices in the boardValid array for later.
	G3VectorFrameObjectConstPtr board_id;
	board_id = boost::dynamic_pointer_cast<G3VectorFrameObject>(
	    bolos->at("board_id"));
	for (size_t i = 0; i < board_id->size(); i++) {
		int board_num = atoi(boost::dynamic_pointer_cast<G3String>
                    ((*board_id)[i])->value.c_str());
		board_id_map_[board_num] = i;
	}

	// Switch on old vs. new data, extracting wiring information from the
	// pipeline iff it is absent in the archive files.
	if (bolos->find("dfml_addr") == bolos->end()) {
		if (!cached_wiring_map_)
			log_fatal("No wiring information (dfml_addr registers) "
			    "in the data and no wiring map to make up for it. "
			    "Please insert a wiring map into the pipeline "
			    "ahead of the first GcpSlow frame containing the "
			    "detector ID->hardware path information.");

		for (auto i : *phys_id_reg) {
			std::string id =
			    boost::dynamic_pointer_cast<G3String>(i)->value;
			if (cached_wiring_map_->find(id) ==
			    cached_wiring_map_->end()) {
				// Add a bogus channel as a sentinel for
				// channels not in the wiring map.
				DfMuxChannelMapping channel_info;
				channel_info.board_serial = -2;
				channels_.push_back(channel_info);
			} else {
				const DfMuxChannelMapping &channel_info =
				    cached_wiring_map_->at(id);
				channels_.push_back(channel_info);
				if (channel_info.channel > max_channel_count_)
					max_channel_count_ =
					    channel_info.channel + 1;
			}
		} 
		return G3FramePtr();
	} else {
		if (cached_wiring_map_)
			log_warn("External wiring frame found in data stream, "
			    "but archive data contains wiring information. "
			    "Emitting a new wiring frame with the archive "
			    "information.");
	}

	dfml_addr = boost::dynamic_pointer_cast<G3VectorFrameObject>(
	    bolos->at("dfml_addr"));
	mbi_addr = boost::dynamic_pointer_cast<G3VectorFrameObject>(
	    bolos->at("mbi_addr"));

	g3_assert(dfml_addr->size() == mbi_addr->size());
	g3_assert(dfml_addr->size() == phys_id_reg->size());


	for (size_t i = 0; i < dfml_addr->size(); i++) {
		DfMuxChannelMapping channel_info;
		std::string channel_path, phys_id, logical_id;
		std::vector<int> channel_parts;

		// Fill fields that don't apply for legacy mux
		channel_info.crate_serial = -1;
		channel_info.board_slot = -1;

		// Parse channel ID (board_module_channel)
		channel_path = boost::dynamic_pointer_cast<G3String>
                    ((*mbi_addr)[i])->value;

		// Special-case the calibrator board
		if (channel_path == "CALIB") {
			// Flag this for later to keep indexing
			channel_info.board_serial = -1;
			channel_info.module = 1;
			channel_info.channel = 0;
			channel_info.board_ip = FAKE_BOARD_IP(255);
			channels_.push_back(channel_info);
			continue;
		}

		boost::char_separator<char> sep("_");
		boost::tokenizer<boost::char_separator<char> > tok(
		    channel_path, sep);
		for (auto j = tok.begin(); j != tok.end(); j++)
			channel_parts.push_back(atoi(j->c_str()));

		// Assign to structure. GCP channel IDs are 1-indexed. 3G's are
		// 0-indexed, so shift by 1.
		g3_assert(channel_parts.size() == 3);
		channel_info.board_serial = channel_parts[0];
		channel_info.module = channel_parts[1] - 1;
		channel_info.channel = channel_parts[2] - 1;
		g3_assert(channel_info.channel >= 0);

		// Guess an IP address. This isn't important.
		channel_info.board_ip = FAKE_BOARD_IP(channel_parts[0]);

		// Get channel ID strings
		phys_id = boost::dynamic_pointer_cast<G3String>
                    ((*phys_id_reg)[i])->value;
		logical_id = boost::dynamic_pointer_cast<G3String>
                    ((*dfml_addr)[i])->value;
		(*wiring)[logical_id] = channel_info;
		(*phys_ids)[logical_id] = phys_id;
		channels_.push_back(channel_info);

		// Book-keeping assistance: we need to know how many channels
		// these boards have (it is 16, but pretend we don't know for
		// robustness) in order to build DfMuxSample objects later.
		if (channel_info.channel >= max_channel_count_)
			max_channel_count_ = channel_info.channel + 1;
	}

	// Now into the frame
	G3FramePtr wiringframe(new G3Frame(G3Frame::Wiring));
	wiringframe->Put("WiringMap", wiring);
	wiringframe->Put("PhysicalBoloIDs", phys_ids);
	wiringframe->Put("ReadoutSystem", G3StringPtr(new G3String("DfMux")));

	log_debug("%zd channels entered into wiring map", wiring->size());
	log_debug("Maximum channel count seen is %zd", max_channel_count_);

	return wiringframe;
}

void
GCPMuxDataDecoder::Process(G3FramePtr frame, std::deque<G3FramePtr> &out_queue)
{
	// For ancient SPT data without stored wiring information, we can
	// use an external wiring map, which we need to cache.
	if (frame->type == G3Frame::Wiring && frame->Has("WiringMap"))
		cached_wiring_map_ = frame->Get<DfMuxWiringMap>("WiringMap");

	// Otherwise we only operate on GCP slow frames
	if (frame->type != G3Frame::GcpSlow) {
		out_queue.push_back(frame);
		return;
	}

	if (!hwm_emitted_) {
		G3FramePtr wf = EmitWiringMap(frame);
		if (wf)
			out_queue.push_back(wf);
		hwm_emitted_ = true;
	}

	out_queue.push_back(frame);

	// Now build a sequence of Timepoint frames
	G3MapFrameObjectConstPtr recv, bolos;
	G3VectorFrameObjectConstPtr times, samplesValid_wrap, adcI, adcQ;
	G3VectorFrameObjectConstPtr boardValid;
	G3VectorIntConstPtr samplesValid;
	std::vector<G3VectorIntConstPtr> Icounts, Qcounts, valid_flag;

	recv = frame->Get<G3MapFrameObject>("receiver");
	bolos = boost::dynamic_pointer_cast<G3MapFrameObject>(
	    recv->at("bolometers"));

	times = boost::dynamic_pointer_cast<G3VectorFrameObject>(
	    bolos->at("utc"));
	times = boost::dynamic_pointer_cast<G3VectorFrameObject>(
	    times->at(0));
	samplesValid_wrap = boost::dynamic_pointer_cast<G3VectorFrameObject>(
	    bolos->at("samplesValid"));
	samplesValid = boost::dynamic_pointer_cast<G3VectorInt>(
	    samplesValid_wrap->at(0));
	boardValid = boost::dynamic_pointer_cast<G3VectorFrameObject>(
	    bolos->at("boardValid"));
	adcI = boost::dynamic_pointer_cast<G3VectorFrameObject>(
	    bolos->at("adcCountsI"));
	adcQ = boost::dynamic_pointer_cast<G3VectorFrameObject>(
	    bolos->at("adcCountsQ"));

	g3_assert(adcI->size() == adcQ->size());
	g3_assert(adcI->size() == channels_.size());

	// Reorganize data and do casts to make this easier to access
	for (size_t i = 0; i < adcI->size(); i++) {
		G3VectorIntConstPtr I = boost::dynamic_pointer_cast<G3VectorInt>
		    (adcI->at(i));
		G3VectorIntConstPtr Q = boost::dynamic_pointer_cast<G3VectorInt>
		    (adcQ->at(i));
		Icounts.push_back(I);
		Qcounts.push_back(Q);
	}

	g3_assert(boardValid->size() == board_id_map_.size());

	for (size_t i = 0; i < boardValid->size(); i++) {
		G3VectorIntConstPtr V = boost::dynamic_pointer_cast<G3VectorInt>
		    (boardValid->at(i));
		valid_flag.push_back(V);
	}

	// Now make some frames
	for (size_t i = 0; i < samplesValid->size(); i++) {
		if (!samplesValid->at(i))
			continue;

		G3FramePtr out(new G3Frame(G3Frame::Timepoint));

		out->Put("EventHeader", times->at(i));

		DfMuxMetaSamplePtr metasample(new DfMuxMetaSample);
		G3DoublePtr calibrator(new G3Double(NAN));
		G3TimeStamp timestamp = boost::dynamic_pointer_cast<G3Time>(
		    times->at(i))->time;

		for (size_t j = 0; j < channels_.size(); j++) {
			const DfMuxChannelMapping &chan = channels_[j];

			// Special-case the calibrator timestream
			if (chan.board_serial == -1) {
				*calibrator = (*Icounts[j])[i] >> 8;
				continue;
			}

			// Skip virtual channels
			if (chan.board_serial == -2)
				continue;

			DfMuxSamplePtr sample =
			    (*metasample)[chan.board_serial][chan.module];

			// Lazily instantiate. Submaps are lazily made by the
			// [][] above (std::map operator [] creates elements
			// if they do not already exist).
			if (!sample) {
				sample = DfMuxSamplePtr(new DfMuxSample(
				    timestamp, 2*max_channel_count_));
				(*metasample)[chan.board_serial][chan.module] =
				    sample;
			}

			// GCP aligns the 24-bit samples to the left. 3G
			// aligns to the right. Use the 3G convention.
			(*sample)[chan.channel*2] = (*Icounts[j])[i] / 256;
			(*sample)[chan.channel*2+1] = (*Qcounts[j])[i] / 256;
		}

		// Set nmodules to avoid triggering warnings
		for (auto j = metasample->begin(); j != metasample->end(); j++)
			j->second.nmodules = j->second.size();

		// Process boardValid: delete boards that are invalid. This
		// corresponds to the standard missing board signalling in
		// 3G.
		for (auto j = board_id_map_.begin(); j != board_id_map_.end();
		    j++) {
			if (!(*valid_flag[j->second])[i])
				metasample->erase(j->first);
		}

		out->Put("DfMux", metasample);
		out->Put("Calibrator", calibrator);

		out_queue.push_back(out);
	}
}

