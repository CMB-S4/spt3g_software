#include <pybindings.h>

#include <dfmux/DfMuxBuilder.h>
#include <dfmux/HardwareMap.h>

#include <math.h>
#include <string>

#include <netcdf.h>

class NetCDFDump : public G3Module {
public:
	NetCDFDump(std::string filename);
	virtual ~NetCDFDump();
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
private:
	int nc_file_id_, nc_dim_;
	int nc_time_id_;
	int nc_current_sample_;
	DfMuxWiringMapConstPtr hwm_;
	std::map<std::string, std::pair<int, int> > nc_ids_;
	
	SET_LOGGER("NetCDFDump");
};

NetCDFDump::NetCDFDump(std::string filename)
{
	int err;

	err = nc_create(filename.c_str(), NC_64BIT_OFFSET | NC_SHARE,
	    &nc_file_id_);
	if (err)
		log_fatal("Error opening %s for writing: %s", filename.c_str(),
		    nc_strerror(err));
	nc_def_dim(nc_file_id_, "time", NC_UNLIMITED, &nc_dim_);
	nc_set_fill(nc_file_id_, NC_NOFILL, NULL);
	nc_current_sample_ = 0;

	nc_def_var(nc_file_id_, "Time", NC_DOUBLE, 1, &nc_dim_, &nc_time_id_);
}

NetCDFDump::~NetCDFDump()
{
	nc_close(nc_file_id_);
}

void NetCDFDump::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
{
	out.push_back(frame);

	if (frame->type == G3Frame::Wiring) {
		hwm_ = frame->Get<DfMuxWiringMap>("WiringMap");

		// Sort channel list alphabetically to make things easier in KST
		std::vector<std::string> chanlist;
		for (auto chan = hwm_->begin(); chan != hwm_->end(); chan++)
			chanlist.push_back(chan->first);
		std::sort(chanlist.begin(), chanlist.end());

		for (auto chan = chanlist.begin(); chan != chanlist.end(); chan++) {
			std::string boloname = *chan;
			std::pair<int, int> &ids = nc_ids_[boloname];
			int err;

			// NetCDF doesn't allow / in names; replace with _
			std::replace(boloname.begin(), boloname.end(),
			    '/', '_');

			err = nc_def_var(nc_file_id_, (boloname + "_I").c_str(),
			    NC_INT, 1, &nc_dim_, &ids.first);
			if (err != NC_NOERR)
				log_fatal("Error creating I column for "
				    "detector %s: %s", boloname.c_str(),
				    nc_strerror(err));

			err = nc_def_var(nc_file_id_, (boloname + "_Q").c_str(),
			    NC_INT, 1, &nc_dim_, &ids.second);
			if (err != NC_NOERR)
				log_fatal("Error creating Q column for "
				    "detector %s: %s", boloname.c_str(),
				    nc_strerror(err));
		}

		nc_enddef(nc_file_id_);
	}

	if (frame->type == G3Frame::Timepoint) {
		DfMuxMetaSampleConstPtr metasamp =
		    frame->Get<DfMuxMetaSample>("DfMux");
		size_t count = 1;
		size_t start = nc_current_sample_;
		int32_t i, q;

		if (frame->Has<G3Time>("EventHeader")) {
			double timestamp =
			    frame->Get<G3Time>("EventHeader")->time/G3Units::s;

			nc_put_vara(nc_file_id_, nc_time_id_, &start, &count,
			    &timestamp);
		}
			

		for (auto chan = hwm_->begin(); chan != hwm_->end(); chan++) {
			const std::pair<int, int> &ids = nc_ids_[chan->first];
			auto board = metasamp->find(chan->second.board_serial);
			if (board == metasamp->end())
				continue;

			auto module = board->second.find(chan->second.module);
			if (module == board->second.end())
				continue;

			if (module->second->size()/2 <=
			    size_t(chan->second.channel)) {
				log_fatal("Board %d, module %d only has %zd "
				    "channels, but trying to read %d",
				    chan->second.board_serial,
				    chan->second.module,
				    module->second->size()/2,
				    chan->second.channel);
			}

			i = (*module->second)[chan->second.channel*2];
			q = (*module->second)[chan->second.channel*2 + 1];

			nc_put_vara(nc_file_id_, ids.first, &start, &count, &i);
			nc_put_vara(nc_file_id_, ids.second, &start, &count, &q);
		}

		nc_current_sample_++;
	}
}

EXPORT_G3MODULE("dfmux", NetCDFDump, init<std::string>(args("filename")),
    "Writes DfMux streamer data to a NetCDF file");

