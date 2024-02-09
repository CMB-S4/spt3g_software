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
#include <calibration/BoloProperties.h>

class HitsBinner : public G3Module {
public:
	HitsBinner(std::string output_map_id, const G3SkyMap &stub_map,
	    std::string pointing, std::string timestreams,
	    std::string bolo_properties_name, boost::python::object map_per_scan);
	virtual ~HitsBinner() {}

	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);

private:
	void BinHits(const BolometerProperties &bp, const G3VectorQuat &pointing,
	    G3SkyMapPtr H);

	std::string output_id_;
	std::string pointing_;
	std::string timestreams_;
	std::string boloprops_name_;
	int map_per_scan_;
	boost::python::object map_per_scan_callback_;

	G3SkyMapPtr H_;
	G3Time start_, stop_;

	BolometerPropertiesMapConstPtr boloprops_;

	SET_LOGGER("HitsBinner");
};

EXPORT_G3MODULE("maps", HitsBinner,
    (init<std::string, const G3SkyMap &, std::string, std::string, std::string,
     object>((arg("map_id"), arg("stub_map"), arg("pointing"),
     arg("timestreams"), arg("bolo_properties_name")="BolometerProperties",
     arg("map_per_scan")=false))),
"HitsBinner(map_id, stub_map, pointing, timestreams, bolo_properties_name=\"BolometerProperties\", map_per_scan=False)\n"
"\n"
"Bins up pointing for a set of channels into a hits map with properties (projection, etc.) \n"
"specified by the stub_map argument. \n\n"
"Parameters\n"
"----------\n"
"map_id : string\n"
"    Will be set as the \"Id\" element of the output frame from this module.\n"
"stub_map : G3SkyMap\n"
"    Template of the map in which to accumulate timestream data. All \n"
"    parameters of the output map (projection, boundaries, pixel size, etc.) \n"
"    are copied from this map, which is not modified.\n"
"pointing : string (frame object G3TimestreamQuat)\n"
"    Name of a frame object containing the boresight pointing quaternion \n"
"    timestream. Must be present in every Scan frame.\n"
"timestreams : string (frame object G3TimestreamMap)\n"
"    Name of a frame object containing the timestreams to be binned into the \n"
"    output map. Must exist in every Scan frame that is to be binned. \n"
"    This key is used only to determine the subset of channels to bin into \n"
"    the hits map.\n"
"bolo_properties_name : string, optional (frame object BolometerPropertiesMap)\n"
"    Name of a bolometer properties map containing detector pointing offsets\n"
"    from boresight. These are usually named \"BolometerProperties\", the \n"
"    default, and this parameter only need be set if you wish to use \n"
"    alternative values.\n"
"map_per_scan: boolean or callable, optional\n"
"    Defaults to False. If set to True, will make an output map for every \n"
"    input Scan frame, which can be useful for fine-grained time-domain \n"
"    maps. Can also be set to a Python callable for complex situations. In \n"
"    this case, the callable will be called on each Scan frame and passed \n"
"    the frame as an argument. If it returns True, the map binner will emit \n"
"    a map frame with the data since the last emitted map frame but *before* \n"
"    the current Scan; if False, the binner will continue accumulating data.\n"
"\n"
"Emits\n"
"-----\n"
"Emits a frame of type Map at the end of processing containing the output \n"
"map:\n\n"
"Frame (Map) [\n"
"\"Id\" (spt3g.core.G3String) => \"ra0hdec-57.5-150GHz\"\n"
"\"H\" (spt3g.maps.FlatSkyMap) => 2700 x 1500 (45 x 25 deg) ZEA centered at (1350, 749.5) = (0, -57.5 deg) in equatorial coordinates  (Tcmb)\n"
"]\n"
"\n"
"See Also\n"
"--------\n"
"MapBinner, SingleDetectorMapBinner, SingleDetectorBoresightBinner :\n"
"    Alternative map binners.\n"
"FlatSkyMap, HealpixSkyMap :\n"
"    Possible output types.\n"
"\n"
"Examples\n"
"--------\n"
".. code-block:: python\n"
"\n"
"    pipe.Add(\n"
"        maps.HitsBinner,\n"
"        map_id=\"150GHz\",\n"
"        stub_map=map_params,\n"
"        timestreams=\"DeflaggedTimestreams150GHz\",\n"
"        pointing=\"OfflineRaDecRotation\",\n"
"    )\n"
);


HitsBinner::HitsBinner(std::string output_map_id, const G3SkyMap &stub_map,
    std::string pointing, std::string timestreams, std::string bolo_properties_name,
    boost::python::object map_per_scan) :
  output_id_(output_map_id), pointing_(pointing), timestreams_(timestreams),
  boloprops_name_(bolo_properties_name)
{
	H_ = stub_map.Clone(false);
	H_->pol_type = G3SkyMap::None;
	H_->SetPolConv(G3SkyMap::ConvNone);
	H_->units = G3Timestream::None;
	H_->weighted = false;

	if (PyCallable_Check(map_per_scan.ptr())) {
		map_per_scan_callback_ = map_per_scan;
		map_per_scan_ = -1;
	} else {
		map_per_scan_ = boost::python::extract<bool>(map_per_scan)();
	}
}

void
HitsBinner::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
{
	if (frame->Has(boloprops_name_))
		boloprops_ = frame->Get<BolometerPropertiesMap>(
		    boloprops_name_);

	bool emit_map_now = false;
	if (frame->type == G3Frame::EndProcessing)
		emit_map_now = true;
	else if (frame->type == G3Frame::Scan) {
		if (start_.time != 0) { // Data exist
			if (map_per_scan_ >= 0)
				emit_map_now = map_per_scan_;
			else
				emit_map_now = map_per_scan_callback_(frame);
		}
	}
	
	if (emit_map_now) {
		if (start_.time == 0)
			log_error("No valid scan frames found for map %s",
			    output_id_.c_str());

		G3FramePtr out_frame(new G3Frame(G3Frame::Map));
		out_frame->Put("Id", G3StringPtr(new G3String(output_id_)));
		out_frame->Put("H", boost::dynamic_pointer_cast<G3FrameObject>(H_));
		H_ = H_->Clone(false);

		out_frame->Put("StartTime", G3TimePtr(new G3Time(start_)));
		out_frame->Put("StopTime", G3TimePtr(new G3Time(stop_)));
		start_.time = stop_.time = 0;

		out.push_back(out_frame);
	}

	if (frame->type != G3Frame::Scan) {
		out.push_back(frame);
		return;
	}

	if (!frame->Has(timestreams_)) {
		out.push_back(frame);
		return;
	}

	if (!boloprops_)
		log_fatal("Need bolometer properties before detector data "
		    "can be processed.");

	// NB: First error (wrong type) is regarded as more severe
	// than the pointing simply being absent, since it is a
	// configuration error rather than a data error, so it is
	// guaranteed nothing will ever work.
	if (!frame->Has<G3TimestreamQuat>(pointing_) &&
	    frame->Has<G3VectorQuat>(pointing_))
		log_fatal("Pointing %s is a G3VectorQuat, but must "
		    "contain timing information for simulation. "
		    "Please turn it into a G3TimestreamQuat before "
		    "running this module, for example by adding the "
		    "shim module maps.AddTimingToPointingQuats.",
		    pointing_.c_str());

	G3TimestreamQuatConstPtr pointing =
	    frame->Get<G3TimestreamQuat>(pointing_, false);
	if (!pointing) {
		log_error("Missing pointing %s", pointing_.c_str());
		out.push_back(frame);
		return;
	}

	G3TimestreamMapConstPtr timestreams =
	    frame->Get<G3TimestreamMap>(timestreams_, false);
	if (!timestreams) {
		log_error("Missing timestreams %s", timestreams_.c_str());
		out.push_back(frame);
		return;
	}

	// Update start and stop times for map
	if (start_.time == 0 || start_ > pointing->start)
		start_ = pointing->start;
	if (stop_ < pointing->stop)
		stop_ = pointing->stop;

#ifdef OPENMP_FOUND
	if (omp_get_num_threads() > 1) {
		// G3SkyMap is not thread-safe when sparse. Avoid the cost of locks
		// by making one set of sparse maps per thread per scan (which use
		// essentially zero memory) and then co-adding them at the end of
		// the scan.

		// Create a list of detectors to satisfy OpenMP's need for scalar
		// iteration
		std::vector<std::string> dets;
		for (auto i : *timestreams)
			dets.push_back(i.first);

		// Clamp num_threads to prevent memory balloon?
		#pragma omp parallel
		{
			G3SkyMapPtr scanH;

			scanH = H_->Clone(false);

			#pragma omp for
			for (size_t i = 0; i < dets.size(); i++) {
				BinHits(boloprops_->at(dets[i]), *pointing, scanH);
			}

			#pragma omp critical
			{
				// One thread at a time here
				(*H_) += *scanH;
			}
		}

	} else {
#endif
		for (auto ts : *timestreams) {
			BinHits(boloprops_->at(ts.first), *pointing, H_);
		}
#ifdef OPENMP_FOUND
	}
#endif

	out.push_back(frame);
}

void
HitsBinner::BinHits(const BolometerProperties &bp, const G3VectorQuat &pointing,
    G3SkyMapPtr H)
    // Do shared pointers add too much overhead here? Probably not...
{
	// Get per-detector pointing timestream
	auto detpointing = get_detector_pointing_pixels(bp.x_offset, bp.y_offset,
	    pointing, H);

	for (size_t i = 0; i < detpointing.size(); i++) {
		(*H)[detpointing[i]] += 1.0;
	}
}
