#ifndef _G3_ARCWRITER_H
#define _G3_ARCWRITER_H

#include <string>
#include <boost/iostreams/filtering_stream.hpp>

#include <G3Module.h>

class G3Writer : public G3Module {
public:
	G3Writer(std::string filename,
	    std::vector<G3Frame::FrameType> streams={}, bool append=false);
	// Writes to file <filename> all frames with types in <streams>.
	// If <streams> is empty (default), writes all frames.

	virtual ~G3Writer();
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
	void Flush();
private:
	std::string filename_;
	boost::iostreams::filtering_ostream stream_;
	std::vector<G3Frame::FrameType> streams_;

	SET_LOGGER("G3Writer");
};

G3_POINTERS(G3Writer);

#endif
