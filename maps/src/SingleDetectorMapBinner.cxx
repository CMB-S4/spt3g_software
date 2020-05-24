#include <pybindings.h>
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

#include <G3Module.h>
#include <G3Timestream.h>
#include <G3Quat.h>
#include <G3Data.h>
#include <G3Map.h>
#include <maps/G3SkyMap.h>
#include <calibration/BoloProperties.h>

class SingleDetectorBoresightBinner : public G3Module {
public:
	SingleDetectorBoresightBinner(const G3SkyMap &stub_map,
	    std::string pointing, std::string timestreams);
	virtual ~SingleDetectorBoresightBinner() {}

	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);

private:
	std::string pointing_;
	std::string timestreams_;

	G3SkyMapPtr template_;
	std::map<std::string, G3SkyMapPtr> maps_;
	G3SkyMapWeightsPtr map_weights_;

	SET_LOGGER("SingleDetectorBoresightBinner");
};

EXPORT_G3MODULE("maps", SingleDetectorBoresightBinner,
    (init<const G3SkyMap &, std::string, std::string>
     ((arg("stub_map"), arg("pointing"), arg("timestreams")))),
"Makes a simple binned map of the sky, in boresight coordinates, for every "
"detector present in the given <timestreams>. The common pointing is specified ""by the <pointing> argument and the map parameters are given by <stub_map>. "
"When processing ends, this module will emit one map frame per detector. "
"Because the maps are unweighted (weights don't make sense for single-detector "
"maps) and the effective pointing for all detectors is identical, the "
"weight/hit map is identical for all detectors and is stored only in the first "
"output frame. This module is intended for use when making maps for "
"detector-pointing calibration. Because of the single stored weight/hit map, "
"the set of detectors in every scan must be identical.");

SingleDetectorBoresightBinner::SingleDetectorBoresightBinner(
    const G3SkyMap &stub_map, std::string pointing, std::string timestreams) :
  pointing_(pointing), timestreams_(timestreams)
{
	template_ = stub_map.Clone(false);
	template_->pol_type = G3SkyMap::T;
	map_weights_ = G3SkyMapWeightsPtr(
	    new G3SkyMapWeights(template_, false));
}

void
SingleDetectorBoresightBinner::Process(G3FramePtr frame,
    std::deque<G3FramePtr> &out_queue)
{

	if (frame->type == G3Frame::EndProcessing) {
		for (auto i : maps_) {
			G3FramePtr out(new G3Frame(G3Frame::Map));
			out->Put("Id", G3StringPtr(new G3String(i.first)));
			out->Put("T",
			  boost::dynamic_pointer_cast<G3FrameObject>(i.second));

			if (map_weights_) {
				out->Put("Wunpol", map_weights_);
				map_weights_.reset();
			}
			out_queue.push_back(out);
		}
		maps_.clear(); // Don't sit on this memory anymore
		out_queue.push_back(frame);
		return;
	}

	if (frame->type != G3Frame::Scan) {
		out_queue.push_back(frame);
		return;
	}

	G3VectorQuatConstPtr pointing =
	    frame->Get<G3VectorQuat>(pointing_);
	if (!pointing) {
		log_error("Missing pointing %s", pointing_.c_str());
		out_queue.push_back(frame);
		return;
	}

	G3TimestreamMapConstPtr timestreams =
	    frame->Get<G3TimestreamMap>(timestreams_);
	if (!timestreams) {
		log_error("Missing timestreams %s", timestreams_.c_str());
		out_queue.push_back(frame);
		return;
	}

	// Initialize units and maps on the first scan frame (when we have
	// no in-progress maps).
	if (maps_.empty()) {
		template_->units = timestreams->GetUnits();
		for (auto i : *timestreams)
			maps_[i.first] = template_->Clone(false);
	} else {
		if (template_->units != timestreams->GetUnits())
			log_fatal("Timestreams have units that do not match "
			    "earlier timestreams in the pipeline.");
		// XXX Check for lost detectors here?
	}

	// Calculate common pointing
	g3_assert(timestreams->NSamples() == pointing->size());
	std::vector<size_t> pixels = template_->QuatAnglesToPixels(*pointing);
	for (size_t i = 0; i < pixels.size(); i++)
		(*map_weights_->TT)[i] += 1;

#ifdef OPENMP_FOUND
	// Create a list of detectors to satisfy OpenMP's need for scalar
	// iteration
	std::vector<std::string> dets;
	for (auto i : maps_)
		dets.push_back(i.first);

	// Clamp num_threads to prevent memory balloon?
	#pragma omp parallel for
	for (size_t i = 0; i < dets.size(); i++) {
		const std::string &det = dets[i];
		G3SkyMapPtr m = maps_.at(det);
#else
	for (auto i : maps_) {
		const std::string &det = i.first;
		G3SkyMapPtr m = i.second;
#endif
		G3TimestreamConstPtr ts = timestreams->at(det);
	
		g3_assert(ts->size() == pixels.size());
		for (size_t j = 0; j < ts->size(); j++) {
			(*m)[pixels[j]] += (*ts)[j];
		}

	}

	out_queue.push_back(frame);
}

