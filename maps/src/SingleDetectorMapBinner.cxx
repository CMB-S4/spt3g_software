#include <pybindings.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <G3Module.h>
#include <G3Timestream.h>
#include <G3Quat.h>
#include <G3Data.h>
#include <G3Map.h>
#include <maps/G3SkyMap.h>
#include <maps/pointing.h>
#include <calibration/BoloProperties.h>

class SingleDetectorMapBinner : public G3Module {
public:
	SingleDetectorMapBinner(const G3SkyMap &stub_map,
	    std::string pointing, std::string timestreams,
	    std::string bolo_properties_name);
	virtual ~SingleDetectorMapBinner() {}

	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);

private:
	std::string pointing_;
	std::string timestreams_;
	std::string boloprops_name_;
	BolometerPropertiesMapConstPtr boloprops_;

	G3SkyMapPtr template_;
	std::map<std::string,
	    std::pair<G3SkyMapPtr, G3SkyMapWeightsPtr> > maps_;

#ifdef _OPENMP
	std::vector<std::string> dets_;
#endif

	SET_LOGGER("SingleDetectorMapBinner");
};

EXPORT_G3MODULE("maps", SingleDetectorMapBinner,
    (init<const G3SkyMap &, std::string, std::string, std::string>
     ((arg("stub_map"), arg("pointing"), arg("timestreams"),
       arg("bolo_properties_name")="BolometerProperties"))),
"Makes a simple binned map of the sky, in sky coordinates, for every "
"detector present in the given <timestreams>. Boresight pointing is specified "
"by the <pointing> argument, with per-detector offsets from boresight as "
"stored in the given bolometer properties map. The map parameters are copied "
"from <stub_map>. When processing ends, this module will emit one map frame "
"per detector, including per-detector (unpolarized) weights.");

SingleDetectorMapBinner::SingleDetectorMapBinner(
    const G3SkyMap &stub_map, std::string pointing, std::string timestreams,
    std::string bolo_properties_name) :
  pointing_(pointing), timestreams_(timestreams),
  boloprops_name_(bolo_properties_name)
{
	template_ = stub_map.Clone(false);
	template_->pol_type = G3SkyMap::T;
	template_->pol_conv = G3SkyMap::ConvNone;
}

void
SingleDetectorMapBinner::Process(G3FramePtr frame,
    std::deque<G3FramePtr> &out_queue)
{
	if (frame->Has(boloprops_name_))
		boloprops_ = frame->Get<BolometerPropertiesMap>(
		    boloprops_name_);

	if (frame->type == G3Frame::EndProcessing) {
		for (auto i : maps_) {
			G3FramePtr out(new G3Frame(G3Frame::Map));
			out->Put("Id", G3StringPtr(new G3String(i.first)));
			out->Put("T", std::dynamic_pointer_cast<G3FrameObject>
			    (i.second.first));
			out->Put("Wunpol", i.second.second);
			out_queue.push_back(out);
		}
		maps_.clear(); // Don't sit on this memory anymore
#ifdef _OPENMP
		dets_.clear();
#endif
		out_queue.push_back(frame);
		return;
	}

	if (frame->type != G3Frame::Scan) {
		out_queue.push_back(frame);
		return;
	}

	if (!boloprops_)
		log_fatal("Need bolometer properties before detector data "
		    "can be processed.");

	G3VectorQuatConstPtr pointing =
	    frame->Get<G3VectorQuat>(pointing_, false);
	if (!pointing) {
		log_error("Missing pointing %s", pointing_.c_str());
		out_queue.push_back(frame);
		return;
	}

	G3TimestreamMapConstPtr timestreams =
	    frame->Get<G3TimestreamMap>(timestreams_, false);
	if (!timestreams) {
		log_error("Missing timestreams %s", timestreams_.c_str());
		out_queue.push_back(frame);
		return;
	}
	g3_assert(timestreams->NSamples() == pointing->size());

	// Initialize units and maps on the first scan frame (when we have
	// no in-progress maps) or if there might be new detectors.
	if (maps_.empty() || maps_.size() != timestreams->size()) {
		template_->units = timestreams->GetUnits();
		for (auto i : *timestreams) {
			if (maps_.find(i.first) != maps_.end())
				continue;

			maps_[i.first].first = template_->Clone(false);
			maps_[i.first].second = G3SkyMapWeightsPtr(
			    new G3SkyMapWeights(template_));
#ifdef _OPENMP
			// Create a list of detectors to satisfy OpenMP's need
			// for scalar iteration.
			dets_.push_back(i.first);
#endif
		}
	} else {
		if (template_->units != timestreams->GetUnits())
			log_fatal("Timestreams have units that do not match "
			    "earlier timestreams in the pipeline.");
	}


#ifdef _OPENMP
	#pragma omp parallel for
	for (size_t i = 0; i < dets_.size(); i++) {
		const std::string &det = dets_[i];
		auto mw = maps_.at(det);
		G3SkyMapPtr m = mw.first;
		G3SkyMapPtr w = mw.second->TT;
#else
	for (auto i : maps_) {
		const std::string &det = i.first;
		G3SkyMapPtr m = i.second.first;
		G3SkyMapPtr w = i.second.second->TT;
#endif

		G3TimestreamConstPtr ts;
		try {
			ts = timestreams->at(det);
		} catch (std::out_of_range &e) {
			continue;
		}
	
		// Get per-detector pointing timestream
		auto bp = boloprops_->find(det);
		if (bp == boloprops_->end())
			log_fatal("Missing bolometer properties for %s",
			    det.c_str());
		auto pixels = get_detector_pointing_pixels(
		    bp->second.x_offset, bp->second.y_offset,
		    *pointing, m);

		g3_assert(ts->size() == pixels.size());
		for (size_t j = 0; j < ts->size(); j++) {
			(*m)[pixels[j]] += (*ts)[j];
			(*w)[pixels[j]] += 1;
		}

	}

	out_queue.push_back(frame);
}

