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

class MapTODMasker : public G3Module {
public:
	MapTODMasker(std::string pointing, std::string timestreams,
	    G3SkyMapMaskConstPtr mask, std::string tod_mask,
	    std::string bolo_properties_name);
	virtual ~MapTODMasker() {}

	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);

private:
	std::string pointing_;
	std::string timestreams_;
	G3SkyMapMaskConstPtr mask_;
	std::string output_;

	std::string boloprops_name_;
	BolometerPropertiesMapConstPtr boloprops_;

	SET_LOGGER("MapTODMasker");
};

EXPORT_G3MODULE("maps", MapTODMasker,
    (init<std::string, std::string, G3SkyMapMaskConstPtr, std::string,
     std::string>((arg("pointing"), arg("timestreams"), arg("mask"),
     arg("tod_mask")="FilterMask",
     arg("bolo_properties_name")="BolometerProperties"))),
"Creates a filter mask <tod_mask> to indicate to later TOD filtering that "
"a detector has passed over a \"bad\" part of the sky. This is used, for "
"example, to avoid fitting a polynomial to a timestream while the detector is "
"looking at a bright point source. All timestream values for the detectors in "
"<timestreams> are set to be ignored (true in the output <tod_mask>) if they "
"point at pixels of the input <mask> map containing a non-zero value. "
"Boresight pointing is obtained from the quaternion vector specified by "
"<pointing>. Detector pointing offsets and polarization angles and "
"efficiencies are obtained from the specified BolometerPropertiesMap, which "
"can generally be left at its default value.");

MapTODMasker::MapTODMasker(std::string pointing, std::string timestreams,
    G3SkyMapMaskConstPtr mask, std::string tod_mask,
    std::string bolo_properties_name) :
  pointing_(pointing), timestreams_(timestreams), mask_(mask),
  output_(tod_mask), boloprops_name_(bolo_properties_name)
{
}

void
MapTODMasker::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
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
	G3MapVectorBoolPtr output(new G3MapVectorBool);

	// Add blank timestreams for all detectors
	for (auto i : *tsm)
		(*output)[i.first].resize(0);

	// Create a list of detectors to satisfy OpenMP's need for scalar
	// iteration and provide stable sorting for later trimming.
	std::vector<std::string> dets;
	for (auto i : *tsm)
		dets.push_back(i.first);

#ifdef OPENMP_FOUND
	#pragma omp parallel for
#endif
	for (size_t i = 0; i < dets.size(); i++) {
		std::vector<bool> &det = output->at(dets[i]);
		const BolometerProperties &bp = boloprops_->at(dets[i]);

		// Get per-detector pointing timestream
		std::vector<double> alpha, delta;
		get_detector_pointing(bp.x_offset, bp.y_offset, *pointing,
		    mask_->Parent()->coord_ref, alpha, delta);
		auto detpointing = mask_->Parent()->AnglesToPixels(alpha, delta);

		det.resize(pointing->size());
		for (size_t j = 0; j < det.size(); j++)
			det[j] = mask_->at(detpointing[j]);

		// Find out if any elements are set, delete entry if not
		bool is_set = false;

		for (size_t j = 0; j < det.size(); j++) {
			is_set = (is_set || det[j]);
			if (is_set)
				break;
		}

		if (!is_set)
			det.resize(0);
	}

	// Remove all unset timestream masks. Done serially to avoid
	// need for locking in the loop above.
	for (size_t i = 0; i < dets.size(); i++) {
		auto det_iter = output->find(dets[i]);

		if (det_iter->second.size() == 0)
			output->erase(det_iter);
	}
		

	frame->Put(output_, output);

	out.push_back(frame);
}

