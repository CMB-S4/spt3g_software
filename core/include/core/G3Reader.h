#ifndef _G3_ARCREADER_H
#define _G3_ARCREADER_H

#include <string>
#include <vector>
#include <boost/iostreams/filtering_stream.hpp>

#include <G3Module.h>

class G3Reader : public G3Module {
public:
	G3Reader(std::string filename);
	G3Reader(std::vector<std::string> filenames);

	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);

private:
	void StartFile(std::string path);
	bool prefix_file_;
	std::string cur_file_;
	std::deque<std::string> filename_;
	boost::iostreams::filtering_istream stream_;

	SET_LOGGER("G3Reader");
};

G3_POINTERS(G3Reader);

#endif
