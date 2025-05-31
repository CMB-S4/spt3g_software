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

class __attribute__((visibility("hidden"))) MapBinner : public G3Module {
public:
	MapBinner(std::string output_map_id, const G3SkyMap &stub_map,
	    std::string pointing, std::string timestreams,
	    std::string detector_weights, std::string bolo_properties_name,
	    bool store_weight_map, py::object map_per_scan);
	virtual ~MapBinner() {}

	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);

private:
	void BinTimestream(const G3Timestream &det, double weight,
	    const BolometerProperties &bp, const G3VectorQuat &pointing,
	    G3SkyMapPtr T, G3SkyMapPtr Q, G3SkyMapPtr U, G3SkyMapWeightsPtr W);

	std::string output_id_;
	std::string pointing_;
	std::string timestreams_;
	std::string weights_;
	std::string boloprops_name_;
	int map_per_scan_;
	py::object map_per_scan_callback_;

	bool units_set_;

	G3SkyMapPtr T_, Q_, U_;
	G3SkyMapWeightsPtr map_weights_;
	G3Time start_, stop_;

	BolometerPropertiesMapConstPtr boloprops_;

	SET_LOGGER("MapBinner");
};

PYBINDINGS("maps", scope) {
register_g3module<MapBinner>(scope, "MapBinner",
"Bins up timestreams into a map with properties (projection, etc.) specified\n"
"by the stub_map argument. \n\n"
"Parameters\n"
"----------\n"
"map_id : string\n"
"    Will be set as the \"Id\" element of the output frame from this module.\n"
"stub_map : G3SkyMap\n"
"    Template of the map in which to accumulate timestream data. All \n"
"    parameters of the output map (projection, boundaries, pixel size, etc.) \n""    are copied from this map, which is not modified. T-only maps are made \n"
"    if the `pol_conv` property of `stub_map` is None; polarized maps are \n"
"    made if it is set to `IAU` or `Cosmo`.\n"
"pointing : string (frame object G3VectorQuat or G3TimestreamQuat)\n"
"    Name of a frame object containing the boresight pointing quaternion \n"
"    timestream. Must have the same number of elements as the data in \n"
"    `timestreams` and be present in every Scan frame.\n"
"timestreams : string (frame object G3TimestreamMap)\n"
"    Name of a frame object containing the timestreams to be binned into the \n"
"    output map. Must exist in every Scan frame, though may be empty if the \n"
"    frame should be ignored. Units in the output map are taken from the \n"
"    units of the detector timestreams.\n"
"\n"
"detector_weights : string, optional (frame object G3MapDouble)\n"
"    Name of a frame object containing weights to assign data from each \n"
"    detector. This frame object maps detector names (as in `timestreams`)\n"
"    to a scalar weight. If unset or set to an empty string (\"\"), \n"
"    detectors will be weighted equally. If set, every detector in \n"
"    `timestreams` must occur in `detector_weights`, but `detector_weights`\n"
"    may contain additional elements not present in `timestreams`.\n"
"bolo_properties_name : string, optional (frame object BolometerPropertiesMap)\n"
"    Name of a bolometer properties map containing detector pointing offsets\n"
"    from boresight. These are usually named \"BolometerProperties\", the \n"
"    default, and this parameter only need be set if you wish to use \n"
"    alternative values.\n"
"store_weight_map : boolean, optional\n"
"    Defaults to True. If set to False, the output weight map will not be \n"
"    stored, though the output maps will *still be weighted*. It should be \n"
"    set to False if and only if you know what the weights are a priori \n"
"    (e.g. mock observing, where they are the same as the data maps).\n"
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
"frame: G3Frame\n"
"    A frame of type Map at the end of processing containing the output map:\n\n"
"    .. code-block::\n\n"
"        Frame (Map) [\n"
"            \"Id\" (spt3g.core.G3String) => \"ra0hdec-57.5-150GHz\"\n"
"            \"T\" (spt3g.maps.FlatSkyMap) => 2700 x 1500 (45 x 25 deg) ZEA centered \n"
"        at (1350, 749.5) = (0, -57.5 deg) in equatorial coordinates  (Tcmb)\n"
"            \"Q\" (spt3g.maps.FlatSkyMap) => 2700 x 1500 (45 x 25 deg) ZEA centered \n"
"        at (1350, 749.5) = (0, -57.5 deg) in equatorial coordinates  (Tcmb)\n"
"            \"U\" (spt3g.maps.FlatSkyMap) => 2700 x 1500 (45 x 25 deg) ZEA centered \n"
"        at (1350, 749.5) = (0, -57.5 deg) in equatorial coordinates  (Tcmb)\n"
"            \"Wpol\" (spt3g.maps.G3SkyMapWeights) => G3SkyMapWeights\n"
"        ]\n"
"\n"
"    For unpolarized maps, the \"Q\" and \"U\" outputs will be absent and the \n"
"    \"Wpol\" frame object will be named \"Wunpol\".\n"
"\n"
"See Also\n"
"--------\n"
"SingleDetectorMapBinner, SingleDetectorBoresightBinner :\n"
"    Alternative map binners for making single-detector maps.\n"
"MapMockObserver :\n"
"    The inverse of this module. Produces timestreams from an input map.\n"
"FlatSkyMap, HealpixSkyMap :\n"
"    Possible output types.\n"
"\n"
"Examples\n"
"--------\n"
".. code-block:: python\n"
"\n"
"    pipe.Add(\n"
"        maps.MapBinner,\n"
"        map_id=\"150GHz\",\n"
"        stub_map=map_params,\n"
"        timestreams=\"DeflaggedTimestreams150GHz\",\n"
"        pointing=\"OfflineRaDecRotation\",\n"
"        detector_weights=\"TodWeights\",\n"
"    )\n"
)
  .def(py::init<std::string, const G3SkyMap &, std::string, std::string, std::string,
     std::string, bool, py::object>(), py::arg("map_id"), py::arg("stub_map"),
     py::arg("pointing"), py::arg("timestreams"), py::arg("detector_weights")="",
     py::arg("bolo_properties_name")="BolometerProperties",
     py::arg("store_weight_map")=true, py::arg("map_per_scan")=false)
;
};


MapBinner::MapBinner(std::string output_map_id, const G3SkyMap &stub_map,
    std::string pointing, std::string timestreams, std::string detector_weights,
    std::string bolo_properties_name, bool store_weight_map,
    py::object map_per_scan) :
  output_id_(output_map_id), pointing_(pointing), timestreams_(timestreams),
  weights_(detector_weights), boloprops_name_(bolo_properties_name),
  units_set_(false)
{
	T_ = stub_map.Clone(false);
	T_->pol_type = G3SkyMap::T;
	if (store_weight_map)
		map_weights_ = G3SkyMapWeightsPtr(new G3SkyMapWeights(T_));

	if (T_->IsPolarized()) {
		Q_ = stub_map.Clone(false);
		Q_->pol_type = G3SkyMap::Q;
		U_ = stub_map.Clone(false);
		U_->pol_type = G3SkyMap::U;
	}

	if (py::isinstance<py::function>(map_per_scan)) {
		map_per_scan_callback_ = map_per_scan;
		map_per_scan_ = -1;
	} else {
		map_per_scan_ = map_per_scan.cast<bool>();
		map_per_scan_callback_ = py::none();
	}
}

void
MapBinner::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
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
				emit_map_now = map_per_scan_callback_(frame).cast<bool>();
		}
	}
	
	if (emit_map_now) {
		if (start_.time == 0)
			log_error("No valid scan frames found for map %s",
			    output_id_.c_str());

		G3FramePtr out_frame(new G3Frame(G3Frame::Map));
		out_frame->Put("Id", G3StringPtr(new G3String(output_id_)));
		out_frame->Put("T",
		    std::dynamic_pointer_cast<G3FrameObject>(T_));
		T_ = T_->Clone(false);

		if (Q_) {
			out_frame->Put("Q",
			    std::dynamic_pointer_cast<G3FrameObject>(Q_));
			out_frame->Put("U",
			    std::dynamic_pointer_cast<G3FrameObject>(U_));
			Q_ = Q_->Clone(false);
			U_ = U_->Clone(false);
		}

		if (map_weights_) {
			out_frame->Put((Q_) ? "Wpol" : "Wunpol", map_weights_);
			map_weights_ = map_weights_->Clone(false);
		}

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
		log_info("Missing timestreams %s", timestreams_.c_str());
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

	G3TimestreamMapConstPtr timestreams =
	    frame->Get<G3TimestreamMap>(timestreams_, false);
	if (!timestreams) {
		log_error("Missing timestreams %s", timestreams_.c_str());
		out.push_back(frame);
		return;
	}

	G3MapDoubleConstPtr weights;
	if (weights_.size() != 0) {
		weights = frame->Get<G3MapDouble>(weights_, false);
		if (!weights) {
			log_error("Missing weights %s", weights_.c_str());
			out.push_back(frame);
			return;
		}
	}

	// Update start and stop times for map
	if (start_.time == 0 || start_ > timestreams->GetStartTime())
		start_ = timestreams->GetStartTime();
	if (stop_ < timestreams->GetStopTime())
		stop_ = timestreams->GetStopTime();

	if (!units_set_) {
		T_->units = timestreams->GetUnits();
		if (Q_)
			Q_->units = timestreams->GetUnits();
		if (U_)
			U_->units = timestreams->GetUnits();
		units_set_ = true;
	} else {
		if (T_->units != timestreams->GetUnits())
			log_fatal("Timestreams have units that do not match "
			    "earlier timestreams in the pipeline.");
	}

#ifdef _OPENMP
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
			G3SkyMapPtr scanT, scanQ, scanU;
			G3SkyMapWeightsPtr scanW;

			scanT = T_->Clone(false);
			if (Q_)
				scanQ = Q_->Clone(false);
			if (U_)
				scanU = U_->Clone(false);
			if (map_weights_)
				scanW = map_weights_->Clone(false);

			#pragma omp for
			for (size_t i = 0; i < dets.size(); i++) {
				BinTimestream(*timestreams->at(dets[i]),
				    (weights) ? weights->at(dets[i]) : 1,
				    boloprops_->at(dets[i]), *pointing,
				    scanT, scanQ, scanU, scanW);
			}

			#pragma omp critical
			{
				// One thread at a time here
				(*T_) += *scanT;
				if (Q_)
					(*Q_) += *scanQ;
				if (U_)
					(*U_) += *scanU;
				if (map_weights_)
					(*map_weights_) += *scanW;
			}
		}

	} else {
#endif
		for (auto ts : *timestreams) {
			BinTimestream(*ts.second,
			    (weights) ? weights->at(ts.first) : 1,
			    boloprops_->at(ts.first), *pointing,
			    T_, Q_, U_, map_weights_);
		}
#ifdef _OPENMP
	}
#endif

	out.push_back(frame);
}

void
MapBinner::BinTimestream(const G3Timestream &det, double weight,
    const BolometerProperties &bp, const G3VectorQuat &pointing,
    G3SkyMapPtr T, G3SkyMapPtr Q, G3SkyMapPtr U, G3SkyMapWeightsPtr W)
    // Do shared pointers add too much overhead here? Probably not...
{
	// Get per-detector pointing timestream

	auto detpointing = get_detector_pointing_pixels(bp.x_offset, bp.y_offset,
	    pointing, T);

	if (Q) {
		// And polarization coupling
		// XXX: does not handle focal-plane rotation, since it assumes
		// polarization coupling is constant for the whole scan.
		StokesVector pcoupling(bp.pol_angle, bp.pol_efficiency);

		MuellerMatrix mueller;
		if (W) {
			mueller.from_vector(pcoupling);
			mueller *= weight;
		}

		pcoupling *= weight;

		for (size_t i = 0; i < det.size(); i++) {
			(*T)[detpointing[i]] += pcoupling.t*det[i];
			(*Q)[detpointing[i]] += pcoupling.q*det[i];
			(*U)[detpointing[i]] += pcoupling.u*det[i];
			if (W)
				(*W)[detpointing[i]] += mueller;
		}
	} else {
		for (size_t i = 0; i < det.size(); i++) {
			(*T)[detpointing[i]] += weight*det[i];
			if (W)
				(*W->TT)[detpointing[i]] += weight;
		}
	}
}

