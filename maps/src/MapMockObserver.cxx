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

class MapMockObserver : public G3Module {
public:
	MapMockObserver(std::string pointing, std::string timestreams,
	    double band, G3SkyMapConstPtr T, G3SkyMapConstPtr Q,
	    G3SkyMapConstPtr U, std::string bolo_properties_name);
	virtual ~MapMockObserver() {}

	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);

private:
	void BinTimestream(const G3Timestream &det, double weight,
	    const BolometerProperties &bp, const G3VectorQuat &pointing,
	    G3SkyMapPtr T, G3SkyMapPtr Q, G3SkyMapPtr U, G3SkyMapWeightsPtr W);

	std::string pointing_;
	std::string timestreams_;
	double band_;
	G3SkyMapConstPtr T_, Q_, U_;

	std::string boloprops_name_;
	BolometerPropertiesMapConstPtr boloprops_;

	SET_LOGGER("MapMockObserver");
};

EXPORT_G3MODULE("maps", MapMockObserver,
    (init<std::string, std::string, double,
     G3SkyMapConstPtr, G3SkyMapConstPtr, G3SkyMapConstPtr,
     std::string>((arg("pointing"), arg("timestreams"), arg("band"),
     arg("T"), arg("Q")=G3SkyMapConstPtr(), arg("U")=G3SkyMapConstPtr(),
     arg("bolo_properties_name")="BolometerProperties"))),
"Creates a new set of timestreams <timestreams> with the expected mean values "
"on each of the detectors with band <band> of a nearest-neighbor observation "
"of the input maps (T and optionally Q and U). Boresight pointing is obtained "
"from the quaternion vector specified by <pointing>. Detector pointing offsets "
"and polarization angles and efficiencies are obtained from the specified "
"BolometerPropertiesMap, which can generally be left at its default value.");

MapMockObserver::MapMockObserver(std::string pointing, std::string timestreams,
    double band, G3SkyMapConstPtr T, G3SkyMapConstPtr Q,
    G3SkyMapConstPtr U, std::string bolo_properties_name) :
  pointing_(pointing), timestreams_(timestreams), band_(band),
  T_(T), Q_(Q), U_(U), boloprops_name_(bolo_properties_name)
{
	if ((Q_ && !U_) || (U_ && !Q_))
		log_fatal("If simulating polarized maps, pass both Q and U.");
}

// Following utility function copied from old map-maker and needs re-tooling
// to handle boresight rotation.
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
MapMockObserver::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
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

	G3TimestreamQuatConstPtr pointing =
	    frame->Get<G3TimestreamQuat>(pointing_);
	if (!pointing) {
		// NB: First error (wrong type) is regarded as more severe
		// than the pointing simply being absent, since it is a
		// configuration error rather than a data error, so it is
		// guaranteed nothing will ever work.
		if (frame->Has<G3VectorQuat>(pointing_))
			log_fatal("Pointing %s is a G3VectorQuat, but must "
			    "contain timing information for simulation. "
			    "Please turn it into a G3TimestreamQuat before "
			    "running this module, for example by adding the "
			    "shim module maps.AddTimingToPointingQuats.",
			    pointing_.c_str());
		else
			log_error("Missing pointing %s", pointing_.c_str());
		out.push_back(frame);
		return;
	}

	G3TimestreamMapPtr tsm = G3TimestreamMapPtr(new G3TimestreamMap);

	// Add blank timestreams for all detectors
	for (auto i : *boloprops_) {
		if (i.second.band != band_)
			continue;
		(*tsm)[i.first] = G3TimestreamPtr(
		    new G3Timestream(pointing->size()));
	}

#ifdef OPENMP_FOUND
	// Create a list of detectors to satisfy OpenMP's need for scalar
	// iteration
	std::vector<std::string> dets;
	for (auto i : *tsm)
		dets.push_back(i.first);

	#pragma omp parallel for
	for (size_t i = 0; i < dets.size(); i++) {
		G3Timestream &det = *tsm->at(dets[i]);
		const BolometerProperties &bp = boloprops_->at(dets[i]);
#else
	for (auto ts : *tsm) {
		G3Timestream &det = *ts.second;
		const BolometerProperties &bp = boloprops_->at(ts.first);
#endif
		det.units = T_->units;
		det.start = pointing->start;
		det.stop = pointing->stop;

		// Get per-detector pointing timestream
		std::vector<double> alpha, delta;
		get_detector_pointing(bp.x_offset, bp.y_offset, *pointing,
		    T_->coord_ref, alpha, delta);
		auto detpointing = T_->AnglesToPixels(alpha, delta);

		if (Q_) {
			StokesVector pcoupling;
			set_stokes_coupling(bp.pol_angle, bp.pol_efficiency,
			    pcoupling);
			for (size_t i = 0; i < det.size(); i++)
				det[i] =
				    (*T_)[detpointing[i]]*pcoupling.t + 
				    (*Q_)[detpointing[i]]*pcoupling.q +
				    (*U_)[detpointing[i]]*pcoupling.u;
		} else {
			for (size_t i = 0; i < det.size(); i++)
				det[i] = (*T_)[detpointing[i]];
		}
	}

	frame->Put(timestreams_, tsm);

	out.push_back(frame);
}

