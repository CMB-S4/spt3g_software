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

class MapMockObserver : public G3Module {
public:
	MapMockObserver(std::string pointing, std::string timestreams,
	    double band, G3SkyMapConstPtr T, G3SkyMapConstPtr Q,
	    G3SkyMapConstPtr U, std::string bolo_properties_name,
	    bool interp, bool error_on_zero);
	virtual ~MapMockObserver() {}

	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);

private:
	std::string pointing_;
	std::string timestreams_;
	double band_;
	G3SkyMapConstPtr T_, Q_, U_;

	std::string boloprops_name_;
	BolometerPropertiesMapConstPtr boloprops_;

	bool interp_;
	bool error_on_zero_;
	int pol_sign_;

	SET_LOGGER("MapMockObserver");
};

EXPORT_G3MODULE("maps", MapMockObserver,
    (init<std::string, std::string, double,
     G3SkyMapConstPtr, G3SkyMapConstPtr, G3SkyMapConstPtr,
     std::string, bool, bool>((arg("pointing"), arg("timestreams"), arg("band"),
     arg("T"), arg("Q")=G3SkyMapConstPtr(), arg("U")=G3SkyMapConstPtr(),
     arg("bolo_properties_name")="BolometerProperties", arg("interp")=false,
     arg("error_on_zero")=true))),
"MapMockObserver(pointing, timestreams, band, T, Q=None, U=None, bolo_properties_name=\"BolometerProperties\", interp=False)\n"
"\n"
"Creates a new set of timestreams by sampling from an input map.\n\n"
"Parameters\n"
"----------\n"
"pointing : string (frame object G3TimestreamQuat)\n"
"    Name of a frame object containing the boresight pointing quaternion \n"
"    timestream.  Used to construct detector pointing with which to sample \n"
"    the input map(s).  Must be present in every Scan frame.\n"
"timestreams : string (frame object G3TimestreamMap)\n"
"    Name of a frame object that will contain the output timestreams sampled \n"
"    from the input map(s).  Units of the timestreams are taken from the \n"
"    input map.  Will be inserted into every Scan frame.\n"
"band : float\n"
"    Frequency band for which the input map was constructed.  Only detectors \n"
"    whose band matches this value will be included in the output timestreams.\n"
"T : frame object G3SkyMap\n"
"    The input temperature map to be sampled.\n"
"Q, U : frame object G3SkyMap, optional\n"
"    Optional Q and U polarization maps to include in the output timestreams, \n"
"    using the polarization angle and efficiency to determine the appropriate \n"
"    coupling to each Stokes component.\n"
"bolo_properties_name : string, optional (frame object BolometerPropertiesMap)\n"
"    Name of a bolometer properties map containing detector pointing offsets\n"
"    from boresight. These are usually named \"BolometerProperties\", the \n"
"    default, and this parameter only need be set if you wish to use \n"
"    alternative values.\n"
"interp : bool, optional\n"
"    If True, use bilinear interpolation to sample the input map(s).  If \n"
"    False (default), use the nearest-neighbor pixel value.\n"
"error_on_zero : bool, optional\n"
"    If True (default), complain loudly if the simulation scans across zeroes\n"
"    in the input map or out of the input map boundaries. If False, allow\n"
"    this to happen without comment.\n"
"\n"
"See Also\n"
"--------\n"
"MapBinner :\n"
"    The inverse of this module.  Bins timestreams into an output map.\n"
"FlatSkyMap, HealpixSkyMap :\n"
"    Possible input map types\n"
"\n"
"Examples\n"
"--------\n"
".. code-block:: python\n"
"\n"
"    pipe.Add(\n"
"        maps.MapMockObserver,\n"
"        pointing=\"OfflineRaDecRotation\",\n"
"        timestreams=\"CalTimestreams150GHz\",\n"
"        band=150 * core.G3Units.GHz,\n"
"        T=map_frame[\"T\"],\n"
"    )\n"
);

MapMockObserver::MapMockObserver(std::string pointing, std::string timestreams,
    double band, G3SkyMapConstPtr T, G3SkyMapConstPtr Q, G3SkyMapConstPtr U,
    std::string bolo_properties_name, bool interp, bool error_on_zero) :
  pointing_(pointing), timestreams_(timestreams), band_(band),
  T_(T), Q_(Q), U_(U), boloprops_name_(bolo_properties_name),
  interp_(interp), error_on_zero_(error_on_zero)
{
	if ((Q_ && !U_) || (U_ && !Q_))
		log_fatal("If simulating polarized maps, pass both Q and U.");

	if (U_) {
		if (!(U_->IsPolarized()))
			log_fatal("Missing pol_conv");
		pol_sign_ = (U_->pol_conv == G3SkyMap::COSMO) ? -1 : 1;
	}
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

	G3TimestreamMapPtr tsm = G3TimestreamMapPtr(new G3TimestreamMap);

	// Add blank timestreams for all detectors
	for (auto i : *boloprops_) {
		if (i.second.band != band_)
			continue;
		// skip channels with missing offsets
		if (i.second.x_offset != i.second.x_offset)
			continue;
		(*tsm)[i.first] = G3TimestreamPtr(
		    new G3Timestream(pointing->size()));
	}

#ifdef _OPENMP
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
		auto detquats = get_detector_pointing_quats(bp.x_offset, bp.y_offset,
		    *pointing, T_->coord_ref);
		std::vector<size_t> detpointing;
		if (!interp_)
			detpointing = T_->QuatsToPixels(detquats);

		if (Q_) {
			// needs re-tooling to handle boresight rotation.
			StokesVector pcoupling(bp.pol_angle * pol_sign_, bp.pol_efficiency);
			for (size_t i = 0; i < det.size(); i++) {
				if (interp_) {
					std::vector<size_t> pixels;
					std::vector<double> weights;
					T_->GetInterpPixelsWeights(detquats[i], pixels, weights);
					det[i] =
					    T_->GetInterpPrecalc(pixels, weights) * pcoupling.t +
					    Q_->GetInterpPrecalc(pixels, weights) * pcoupling.q +
					    U_->GetInterpPrecalc(pixels, weights) * pcoupling.u;
				} else {
					det[i] =
					    T_->at(detpointing[i]) * pcoupling.t +
					    Q_->at(detpointing[i]) * pcoupling.q +
					    U_->at(detpointing[i]) * pcoupling.u;
				}
			}
		} else {
			for (size_t i = 0; i < det.size(); i++) {
				if (interp_)
					det[i] = T_->GetInterpValue(detquats[i]);
				else
					det[i] = T_->at(detpointing[i]);
			}
		}
		if (error_on_zero_) {
			for (size_t i = 0; i < det.size(); i++) {
				if (det[i] == 0) {
					log_fatal("Scanning across zero-valued "
					    "pixels in input map. Input map is "
					    "likely too small. If this is "
					    "intended, unset parameter "
					    "error_on_zero.");
				}
			}
		}
	}

	frame->Put(timestreams_, tsm);

	out.push_back(frame);
}

