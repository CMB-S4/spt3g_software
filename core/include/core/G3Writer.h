#ifndef _G3_ARCWRITER_H
#define _G3_ARCWRITER_H

#include <string>
#include <boost/iostreams/filtering_stream.hpp>

#include <G3Module.h>

class G3Writer : public G3Module {
public:
	G3Writer(std::string filename);
	virtual ~G3Writer();
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
private:
	std::string filename_;
	boost::iostreams::filtering_ostream stream_;

	SET_LOGGER("G3Writer");
};

G3_POINTERS(G3Writer);

#endif
