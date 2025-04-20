#include <pybindings.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <G3Module.h>
#include <G3Timestream.h>
#include <G3Quat.h>
#include <G3Data.h>
#include <G3Map.h>
#include <G3Units.h>
#include <maps/G3SkyMap.h>
#include <maps/pointing.h>
#include <calibration/BoloProperties.h>

class MapTODPointing : public G3Module {
public:
	MapTODPointing(std::string pointing, std::string timestreams,
	    G3SkyMapConstPtr stub_map, std::string tod_pointing,
	    std::string bolo_properties_name);
	virtual ~MapTODPointing() {}

	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);

private:
	std::string pointing_;
	std::string timestreams_;

	G3SkyMapConstPtr stub_;

	std::string output_;
	std::string boloprops_name_;

	BolometerPropertiesMapConstPtr boloprops_;

	SET_LOGGER("MapTODPointing");
};

PYBINDINGS("maps", scope) {
register_g3module<MapTODPointing>(scope, "MapTODPointing",
"MapTODPointing(pointing, timestreams, stub_map, tod_pointing, bolo_properties_name=\"BolometerProperties\")\n"
"\n"
"Compute pixel pointing for timestreams for a map with properties (projection, etc.) specified\n"
"by the stub_map argument. \n\n"
"The outputs of this module are identical to internally-calculated quantities in the map maker\n"
"and can be used for debugging map-making or studying individual-detector pixel coverage.\n"
"The outputs are not consumed by any stage of the map-making pipeline.\n\n"
"Parameters\n"
"----------\n"
"pointing : string (frame object G3VectorQuat or G3TimestreamQuat)\n"
"    Name of a frame object containing the boresight pointing quaternion \n"
"    timestream. Must have the same number of elements as the data in \n"
"    `timestreams` and be present in every Scan frame.\n"
"timestreams : string (frame object G3TimestreamMap)\n"
"    Name of a frame object that contains the detector timestreams for \n"
"    the channels for which to compute pointing timestreams.\n"
"stub_map : G3SkyMap\n"
"    Template of the map for which to compute pixel pointing.\n"
"tod_pointing : string\n"
"    The key to which the pixel pointing G3MapVectorInt object is\n"
"    stored in the output Scan frames from this module.\n"
"bolo_properties_name : string, optional (frame object BolometerPropertiesMap)\n"
"    Name of a bolometer properties map containing detector pointing offsets\n"
"    from boresight. These are usually named \"BolometerProperties\", the \n"
"    default, and this parameter only need be set if you wish to use \n"
"    alternative values.\n"
"\n"
"Emits\n"
"-----\n"
"Adds a new key to every Scan frame containing the output pointing data:\n\n"
"Frame (Scan) [\n"
"...\n"
"\"OnlineRaDecRotation\" (spt3g.core.G3TimestreamQuat) => 15716 quaternions at 152.6 Hz\n"
"\"RawTimestreams_I\" (spt3g.core.G3TimestreamMap) => Timestreams from 14112 detectors\n"
"\"TodPointing\" (spt3g.core.G3MapVectorInt) => 14112 elements\n"
"...\n"
"]\n"
"\n"
"See Also\n"
"--------\n"
"MapTODMasker :\n"
"    Projection of sky masks to detector timestreams.\n"
"MapBinner, SingleDetectorMapBinner, SingleDetectorBoresightBinner :\n"
"    Map binners for making coadded and single-detector maps.\n"
"MapMockObserver :\n"
"    The inverse of MapBinner. Produces timestreams from an input map.\n"
"FlatSkyMap, HealpixSkyMap :\n"
"    Possible map types.\n"
"\n"
"Examples\n"
"--------\n"
".. code-block:: python\n"
"\n"
"    pipe.Add(\n"
"        maps.MapTODPointing,\n"
"        pointing=\"OnlineRaDecRotation\",\n"
"        timestreams=\"RawTimestreams_I\",\n"
"        stub_map=map_params,\n"
"        tod_pointing=\"TodPointing\",\n"
"    )\n"
)
  .def(py::init<std::string, std::string, G3SkyMapConstPtr, std::string, std::string>((
     py::arg("pointing"), py::arg("timestreams"), py::arg("stub_map"), py::arg("tod_pointing"),
     py::arg("bolo_properties_name")="BolometerProperties")
     ))
;
};


MapTODPointing::MapTODPointing(std::string pointing, std::string timestreams,
    G3SkyMapConstPtr stub_map, std::string tod_pointing,
    std::string bolo_properties_name) :
  pointing_(pointing), timestreams_(timestreams), stub_(stub_map),
  output_(tod_pointing), boloprops_name_(bolo_properties_name)
{
}

void
MapTODPointing::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
{
	if (frame->Has(boloprops_name_))
		boloprops_ = frame->Get<BolometerPropertiesMap>(
		    boloprops_name_);

	if (frame->type != G3Frame::Scan) {
		out.push_back(frame);
		return;
	}

	if (!boloprops_)
		log_fatal("Need bolometer properties before detector data "
		    "can be processed.");

	G3VectorQuatConstPtr pointing =
	    frame->Get<G3VectorQuat>(pointing_, false);
	if (!pointing) {
		log_error("Missing pointing %s", pointing_.c_str());
		out.push_back(frame);
		return;
	}

	G3TimestreamMapConstPtr tsm =
	    frame->Get<G3TimestreamMap>(timestreams_);
	G3MapVectorIntPtr output(new G3MapVectorInt);

	// Add blank timestreams for all detectors and create a list
	// of detectors to satisfy OpenMP's need for scalar iteration
	size_t sz = pointing->size();
	std::vector<std::string> dets;
	for (auto i : *tsm) {
		(*output)[i.first].resize(sz);
		dets.push_back(i.first);
	}

#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (size_t i = 0; i < dets.size(); i++) {
		auto &det = output->at(dets[i]);
		const BolometerProperties &bp = boloprops_->at(dets[i]);

		// Get per-detector pointing timestream
		auto detpointing = get_detector_pointing_pixels(bp.x_offset, bp.y_offset,
		    *pointing, stub_);

		for (size_t j = 0; j < sz; j++)
			det[j] = detpointing[j];
	}

	frame->Put(output_, output);

	out.push_back(frame);
}
