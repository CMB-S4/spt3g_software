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

class MapBinner : public G3Module {
public:
	MapBinner(std::string output_map_id, const G3SkyMap &stub_map,
	    std::string pointing, std::string timestreams,
	    std::string detector_weights, std::string bolo_properties_name,
	    bool store_weight_map);
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

	bool units_set_;

	G3SkyMapPtr T_, Q_, U_;
	G3SkyMapWeightsPtr map_weights_;

	BolometerPropertiesMapConstPtr boloprops_;

	SET_LOGGER("MapBinner");
};

EXPORT_G3MODULE("maps", MapBinner,
    (init<std::string, const G3SkyMap &, std::string, std::string, std::string,
     std::string, bool>((arg("map_id"), arg("stub_map"), arg("pointing"),
     arg("timestreams"), arg("detector_weights")="",
     arg("bolo_properties_name")="BolometerProperties",
     arg("store_weight_map")=true)),
"Bins up timestreams into a map with properties (projection, etc.) specified "
"by the <stub_map> argument. Whether the output map is polarized is determined "
"by the value of the pol_conv property of <stub_map>: T-only maps are made if "
"set to None, T/Q/U are made if set to IAU or Cosmo. If <detector_weights> is "
"set, the given frameobject (a G3MapDouble indexed by detector ID) will be "
"used to weight the detectors when mapped; otherwise, detectors will be "
"weighted equally. Boresight pointing is obtained from the quaternion vector "
"specified by <pointing>. Detector pointing offsets and polarization angles "
"and efficiencies are obtained from the specified BolometerPropertiesMap, "
"which can generally be left at its default value. If <store_weight_map> is "
"set to False (not the default), the output weight map will not be stored,
though the output maps will *still be weighted*. It should be set to False
if and only if you know what the weights are a priori (e.g. mock observing).");

MapBinner::MapBinner(std::string output_map_id, const G3SkyMap &stub_map,
    std::string pointing, std::string timestreams, std::string detector_weights,
    std::string bolo_properties_name, bool store_weight_map) :
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
}

void
MapBinner::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
{
	if (frame->Has(boloprops_name_))
		boloprops_ = frame->Get<BolometerPropertiesMap>(
		    boloprops_name_);

	if (frame->type == G3Frame::EndProcessing) {
		G3FramePtr out_frame(new G3Frame(G3Frame::Map));
		out_frame->Put("Id", G3StringPtr(new G3String(output_id_)));
		out_frame->Put("T",
		    boost::dynamic_pointer_cast<G3FrameObject>(T_));
		if (Q_) {
			out_frame->Put("Q",
			    boost::dynamic_pointer_cast<G3FrameObject>(Q_));
			out_frame->Put("U",
			    boost::dynamic_pointer_cast<G3FrameObject>(U_));
		}

		if (map_weights_)
			out_frame->Put((Q_) ? "Wpol" : "Wunpol", map_weights_);
		out.push_back(out_frame);
		out.push_back(frame);
		return;
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
	G3MapDoubleConstPtr weights = frame->Get<G3MapDouble>(weights_, false);
	if (!timestreams) {
		log_error("Missing timestreams %s", timestreams_.c_str());
		out.push_back(frame);
		return;
	}
	if (!weights) {
		log_error("Missing weights %s", weights_.c_str());
		out.push_back(frame);
		return;
	}

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
#else
	for (auto ts : *timestreams) {
		BinTimestream(*ts.second,
		    (weights) ? weights->at(ts.first) : 1,
		    boloprops_->at(ts.first), *pointing,
		    T_, Q_, U_, map_weights_);
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
		fill_mueller_matrix_from_stokes_coupling(pcoupling, mueller);
		mueller *= weight;
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

