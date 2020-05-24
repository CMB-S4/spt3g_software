#include <pybindings.h>
#include <omp.h>

#include <G3Module.h>
#include <G3Timestream.h>
#include <maps/G3SkyMap.h>
#include <calibration/BoloProperties.h>

class MapBinner : public G3Module {
public:
	MapBinner(std::string output_map_id, const G3SkyMap &stub_map,
	    std::string timestreams, std::string detector_weights,
	    std::string bolo_properties_name);
	virtual ~MapBinner() {}

	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);

private:
	std::string output_id_;
	std::string timestreams_;
	std::string weights_;
	std::string boloprops_name_;

	G3SkyMapPtr T_, Q_, U_;
	G3SkyMapWeightsPtr map_weights_;

	BolometerPropertiesMapConstPtr boloprops_;

	SET_LOGGER("MapBinner");
};

EXPORT_G3MODULE("maps", MapBinner,
    (init<std::string, const G3SkyMap &, std::string, std::string,
     std::string>((arg("map_id"), arg("stub_map"), arg("timestreams"),
     bp::arg("detector_weights")="",
     bp::arg("bolo_properties_name")="BolometerProperties"))),
"Bins up timestreams into a map with properties (projection, etc.) specified "
"by the <stub_map> argument. Whether the output map is polarized is determined "
"by the value of the pol_conv property of <stub_map>: T-only maps are made if "
"set to None, T/Q/U are made if set to IAU or Cosmo. If <detector_weights> is "
"set, the given frameobject (a G3MapDouble indexed by detector ID) will be "
"used to weight the detectors when mapped; otherwise, detectors will be "
"weighted equally. Detector pointing offsets and polarization angles and "
"efficiencies are obtained from the specified BolometerPropertiesMap, which "
"can generally be left at its default value.");

MapBinner::MapBinner(std::string output_map_id, const G3SkyMap &stub_map,
    std::string timestreams, std::string detector_weights,
    std::string bolo_properties_name) :
  output_id_(output_map_id), timestreams_(timestreams),
  weights_(detector_weights), boloprops_name_(bolo_properties_name)
{
	T_ = stub_map.Clone(false);
	T_->pol_type = G3SkyMap::T;
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
		// Emit map frame
		out.push_back(frame);
		return;
	}

	if (frame->type != G3Frame::Scan) {
		out.push_back(frame);
		return;
	}

	G3SkyMapPtr scanT, scanQ, scanU;
	G3SkyMapWeightsPtr scanW;

	struct chunk {
		G3SkyMapPtr scanT, scanQ, scanU;
		G3SkyMapWeightsPtr scanW;
	};
	std::vector<struct chunk> chunk_array;

	G3TimestreamMapConstPtr timestreams =
	    frame->Get<G3TimestreamMap>(timestreams_);
	if (!timestreams) {
		log_warn("Missing timestreams %s", timestreams_.c_str());
		out.push_back(frame);
		return;
	}

	// Create a list of detectors to satisfy OpenMP's need for scalar
	// iteration
	std::vector<std::string> dets;
	for (auto i : *timestreams)
		dets.push_back(i.first);

	// Clamp num_threads to prevent memory balloon?
	#pragma omp parallel private(scanT, scanQ, scanU, scanW) shared(chunk_array)
	{
		#pragma omp single
		if (omp_get_thread_num() == 0)
			chunk_array.resize(omp_get_num_threads());
		#pragma omp barrier

		chunk_array[omp_get_thread_num()].scanT = scanT =
		    T_->Clone(false);
		if (Q_)
			chunk_array[omp_get_thread_num()].scanQ =
			    scanQ = Q_->Clone(false);
		if (U_)
			chunk_array[omp_get_thread_num()].scanU =
			    scanU = U_->Clone(false);
		chunk_array[omp_get_thread_num()].scanW = scanW =
		    map_weights_->Clone(false);

		#pragma omp for
		for (size_t i = 0; i < dets.size(); i++) {
			printf("Detector %zd: %s\n", i, dets[i].c_str());
		}
	}
	for (auto i : chunk_array) {
		(*T_) += *i.scanT;
		if (Q_)
			(*Q_) += *i.scanQ;
		if (U_)
			(*U_) += *i.scanU;
		(*map_weights_) += *i.scanW;
	}

	out.push_back(frame);
}


