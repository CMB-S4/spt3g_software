#include <pybindings.h>
#ifdef OPENMP_FOUND
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

class MapBinner : public G3Module {
public:
	MapBinner(std::string output_map_id, const G3SkyMap &stub_map,
	    std::string pointing, std::string timestreams,
	    std::string detector_weights, std::string bolo_properties_name,
	    bool store_weight_map, boost::python::object map_per_scan);
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
	boost::python::object map_per_scan_callback_;

	bool units_set_;

	G3SkyMapPtr T_, Q_, U_;
	G3SkyMapWeightsPtr map_weights_;
	G3Time start_, stop_;

	BolometerPropertiesMapConstPtr boloprops_;

	SET_LOGGER("MapBinner");
};

EXPORT_G3MODULE("maps", MapBinner,
    (init<std::string, const G3SkyMap &, std::string, std::string, std::string,
     std::string, bool, object>((arg("map_id"), arg("stub_map"),
     arg("pointing"), arg("timestreams"), arg("detector_weights")="",
     arg("bolo_properties_name")="BolometerProperties",
     arg("store_weight_map")=true,arg("map_per_scan")=false))),
"MapBinner(map_id, stub_map, pointing, timestreams, detector_weights, bolo_properties_name=\"BolometerProperties\", store_weight_map=True, map_per_scan=False)\n"
"\n"
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
"Emits a frame of type Map at the end of processing containing the output \n"
"map:\n\n"
"Frame (Map) [\n"
"\"Id\" (spt3g.core.G3String) => \"ra0hdec-57.5-150GHz\"\n"
"\"T\" (spt3g.maps.FlatSkyMap) => 2700 x 1500 (45 x 25 deg) ZEA centered at (1350, 749.5) = (0, -57.5 deg) in equatorial coordinates  (Tcmb)\n"
"\"Q\" (spt3g.maps.FlatSkyMap) => 2700 x 1500 (45 x 25 deg) ZEA centered at (1350, 749.5) = (0, -57.5 deg) in equatorial coordinates  (Tcmb)\n"
"\"U\" (spt3g.maps.FlatSkyMap) => 2700 x 1500 (45 x 25 deg) ZEA centered at (1350, 749.5) = (0, -57.5 deg) in equatorial coordinates  (Tcmb)\n"
"\"Wpol\" (spt3g.maps.G3SkyMapWeights) => G3SkyMapWeights\n"
"]\n"
"\n"
"For unpolarized maps, the \"Q\" and \"U\" outputs will be absent and the \n"
"\"Wpol\" frame object will be named \"Wunpol\".\n"
"\n"
"See Also\n"
"--------\n"
"SingleDetectorMapBinner, SingleDetectorBoresightBinner :\n"
"    Alternative map binners for making single-detector maps.\n"
"MapMockObserver :\n"
"    The inverse of this module. Produces timestreams from aninput map.\n"
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
);


MapBinner::MapBinner(std::string output_map_id, const G3SkyMap &stub_map,
    std::string pointing, std::string timestreams, std::string detector_weights,
    std::string bolo_properties_name, bool store_weight_map,
    boost::python::object map_per_scan) :
  output_id_(output_map_id), pointing_(pointing), timestreams_(timestreams),
  weights_(detector_weights), boloprops_name_(bolo_properties_name),
  units_set_(false)
{
	T_ = stub_map.Clone(false);
	T_->pol_type = G3SkyMap::T;
	if (store_weight_map)
		map_weights_ = G3SkyMapWeightsPtr(new G3SkyMapWeights(T_,
		    T_->GetPolConv() != G3SkyMap::ConvNone));

	if (T_->GetPolConv() != G3SkyMap::ConvNone) {
		Q_ = stub_map.Clone(false);
		Q_->pol_type = G3SkyMap::Q;
		U_ = stub_map.Clone(false);
		U_->pol_type = G3SkyMap::U;
	}

	if (PyCallable_Check(map_per_scan.ptr())) {
		map_per_scan_callback_ = map_per_scan;
		map_per_scan_ = -1;
	} else {
		map_per_scan_ = boost::python::extract<bool>(map_per_scan)();
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
				emit_map_now = map_per_scan_callback_(frame);
		}
	}
	
	if (emit_map_now) {
		G3FramePtr out_frame(new G3Frame(G3Frame::Map));
		out_frame->Put("Id", G3StringPtr(new G3String(output_id_)));
		out_frame->Put("T",
		    boost::dynamic_pointer_cast<G3FrameObject>(T_));
		T_ = T_->Clone(false);

		if (Q_) {
			out_frame->Put("Q",
			    boost::dynamic_pointer_cast<G3FrameObject>(Q_));
			out_frame->Put("U",
			    boost::dynamic_pointer_cast<G3FrameObject>(U_));
			Q_ = Q_->Clone(false);
			U_ = U_->Clone(false);
		}

		if (map_weights_) {
			out_frame->Put((Q_) ? "Wpol" : "Wunpol", map_weights_);
			map_weights_ = G3SkyMapWeightsPtr(new G3SkyMapWeights(
			    T_, T_->GetPolConv() != G3SkyMap::ConvNone));
		}

		out_frame->Put("StartTime", G3TimePtr(new G3Time(start_)));
		out_frame->Put("StopTime", G3TimePtr(new G3Time(stop_)));

		out.push_back(out_frame);
	}

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
	if (start_.time == 0 || start_ > timestreams->start)
		start_ = timestreams->start;
	if (stop_ < timestreams->stop)
		stop_ = timestreams->stop;

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
#ifdef OPENMP_FOUND
	}
#endif

	out.push_back(frame);
}

// Following two utility functions copied from the old map-maker and need to
// be re-tooled -- either moved into the constructors of the relevant objects
// and/or improved to handle things like boresight rotation.
static void
fill_mueller_matrix_from_stokes_coupling(
    const StokesVector &stokes_coupling, MuellerMatrix &pcm)
{
	pcm.tt = stokes_coupling.t * stokes_coupling.t;
	pcm.tq = stokes_coupling.t * stokes_coupling.q;
	pcm.tu = stokes_coupling.t * stokes_coupling.u;
	pcm.qq = stokes_coupling.q * stokes_coupling.q;
	pcm.qu = stokes_coupling.q * stokes_coupling.u;
	pcm.uu = stokes_coupling.u * stokes_coupling.u;
}

static void
set_stokes_coupling(double pol_ang, double pol_eff,
    StokesVector &stokes_coupling)
{
	stokes_coupling.t = 1.0;
	stokes_coupling.q = cos(pol_ang/G3Units::rad *2.)*pol_eff/(2.0-pol_eff);
	stokes_coupling.u = sin(pol_ang/G3Units::rad *2.)*pol_eff/(2.0-pol_eff);

	if (fabs(stokes_coupling.q) < 1e-12)
		stokes_coupling.q = 0;
	if (fabs(stokes_coupling.u) < 1e-12)
		stokes_coupling.u = 0;
}

void
MapBinner::BinTimestream(const G3Timestream &det, double weight,
    const BolometerProperties &bp, const G3VectorQuat &pointing,
    G3SkyMapPtr T, G3SkyMapPtr Q, G3SkyMapPtr U, G3SkyMapWeightsPtr W)
    // Do shared pointers add too much overhead here? Probably not...
{
	// Get per-detector pointing timestream

	std::vector<double> alpha, delta;
	get_detector_pointing(bp.x_offset, bp.y_offset, pointing, T->coord_ref,
	    alpha, delta);

	auto detpointing = T->AnglesToPixels(alpha, delta);

	if (Q) {
		// And polarization coupling
		// XXX: does not handle focal-plane rotation, since it assumes
		// polarization coupling is constant for the whole scan.
		StokesVector pcoupling;
		set_stokes_coupling(bp.pol_angle, bp.pol_efficiency, pcoupling);

		MuellerMatrix mueller;
		if (W) {
			fill_mueller_matrix_from_stokes_coupling(pcoupling, mueller);
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

