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
#include <maps/pointing.h>

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

#ifdef OPENMP_FOUND
	std::vector<std::string> dets_;
#endif

	SET_LOGGER("SingleDetectorBoresightBinner");
};

EXPORT_G3MODULE("maps", SingleDetectorBoresightBinner,
    (init<const G3SkyMap &, std::string, std::string>
     ((arg("stub_map"), arg("pointing"), arg("timestreams")))),
"SingleDetectorBoresightBinner(stub_map, pointing, timestreams)\n"
"\n"
"Makes simple binned maps of the sky, in boresight coordinates, for every \n"
"detector present in the given timestreams. This module is intended for use \n"
"when making maps for detector-pointing calibration.\n"
"\n"
"Parameters\n"
"----------\n"
"stub_map : G3SkyMap\n"
"    Template of the map in which to accumulate timestream data. All \n"
"     parameters of the output map (projection, boundaries, pixel size, etc.)\n"
"     are copied from this map, which is not modified.\n"
"pointing : string (frame object G3VectorQuat or G3TimestreamQuat)\n"
"    Name of a frame object containing the boresight pointing quaternion \n"
"    timestream. Must have the same number of elements as the data in \n"
"    `timestreams` and be present in every Scan frame.\n"
"timestreams : string (frame object G3TimestreamMap)\n"
"    Name of a frame object containing the timestreams to be binned into the \n"
"    output map. Must exist in every Scan frame, though may be empty if the \n"
"    frame should be ignored. Units in the output map are taken from the \n"
"    units of the detector timestreams. Because of the single stored \n"
"    weight/hit map, the set of detectors in every scan must be identical.\n"
"\n"
"Emits\n"
"-----\n"
"At the end of processing, emits one frame of type Map for each detector \n"
"containing the corresponding map. The \"Id\" property of the output frames \n"
"gives the detector ID to which the map corresponds. Because the maps are \n"
"unweighted (weights don't make sense for single-detector maps) and the \n"
"effective pointing for all detectors is identical, the weight/hit map is \n"
"identical for all detectors and is stored only in the first output frame.\n"
"\n"
"Frame (Map) [\n"
"\"Id\" (spt3g.core.G3String) => \"2019.000\"\n"
"\"T\" (spt3g.maps.FlatSkyMap) => 360 x 360 (3 x 3 deg) ZEA centered at (179.5, 179.5) = (0, 0 deg) in local coordinates  (Power)\n"
"\"Wunpol\" (spt3g.maps.G3SkyMapWeights) => G3SkyMapWeights\n"
"]\n"
"Frame (Map) [\n"
"\"Id\" (spt3g.core.G3String) => \"2019.005\"\n"
"\"T\" (spt3g.maps.FlatSkyMap) => 360 x 360 (3 x 3 deg) ZEA centered at (179.5, 179.5) = (0, 0 deg) in local coordinates  (Power)\n"
"]\n"
"\n"
"See Also\n"
"--------\n"
"SingleDetectorMapBinner :\n"
"    Produces single-detector maps, but in sky coordinates, \n"
"    taking into account individual-detector pointing offsets.\n"
"MapBinner :\n"
"    Produces co-added maps from many detectors.\n"
"MapMockObserver :\n"
"    The inverse of MapBinner. Produces timestreams from an input map.\n"
"FlatSkyMap, HealpixSkyMap :\n"
"    Possible output types.\n"
"\n"
"Examples\n"
"--------\n"
".. code-block:: python\n"
"\n"
"    pipe.Add(\n"
"        maps.SingleDetectorBoresightBinner,\n"
"        stub_map=smstub,\n"
"        timestreams='PolyFilteredTimestreams',\n"
"        pointing='OffsetRotation',\n"
"    )\n"
);

SingleDetectorBoresightBinner::SingleDetectorBoresightBinner(
    const G3SkyMap &stub_map, std::string pointing, std::string timestreams) :
  pointing_(pointing), timestreams_(timestreams)
{
	template_ = stub_map.Clone(false);
	template_->pol_type = G3SkyMap::T;
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
#ifdef OPENMP_FOUND
		dets_.clear();
#endif
		out_queue.push_back(frame);
		return;
	}

	if (frame->type != G3Frame::Scan) {
		out_queue.push_back(frame);
		return;
	}

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

	// Initialize units and maps on the first scan frame (when we have
	// no in-progress maps).
	if (maps_.empty()) {
		template_->units = timestreams->GetUnits();
		for (auto i : *timestreams) {
			maps_[i.first] = template_->Clone(false);

			#ifdef OPENMP_FOUND
			// Create a list of detectors to satisfy OpenMP's need
			// for scalar iteration
			dets_.push_back(i.first);
			#endif
		}
		map_weights_ = G3SkyMapWeightsPtr(
		    new G3SkyMapWeights(template_, false));
	} else {
		if (template_->units != timestreams->GetUnits())
			log_fatal("Timestreams have units that do not match "
			    "earlier timestreams in the pipeline.");
		// XXX Check for lost detectors here?
	}

	// Calculate common pointing
	g3_assert(timestreams->NSamples() == pointing->size());
	
	// Conjugate pointing rotation with boresight vector
	auto pointing_vec = get_detector_pointing_quats(0, 0,
	    *pointing, template_->coord_ref);
	auto pixels = template_->QuatsToPixels(pointing_vec);

	for (size_t i = 0; i < pixels.size(); i++)
		(*map_weights_->TT)[pixels[i]] += 1;

#ifdef OPENMP_FOUND
	#pragma omp parallel for
	for (size_t i = 0; i < dets_.size(); i++) {
		const std::string &det = dets_[i];
		G3SkyMapPtr m = maps_.at(det);
#else
	for (auto i : maps_) {
		const std::string &det = i.first;
		G3SkyMapPtr m = i.second;
#endif
		G3TimestreamConstPtr ts;
		try {
			ts = timestreams->at(det);
		} catch (std::out_of_range &e) {
			continue;
		}
	
		g3_assert(ts->size() == pixels.size());
		for (size_t j = 0; j < ts->size(); j++) {
			(*m)[pixels[j]] += (*ts)[j];
		}

	}

	out_queue.push_back(frame);
}

