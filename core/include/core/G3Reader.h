#ifndef _G3_ARCREADER_H
#define _G3_ARCREADER_H

#include <string>
#include <vector>
#include <iostream>

#include <G3Module.h>

class G3Reader : public G3Module {
public:
	G3Reader(const std::string &filename, int n_frames_to_read = -1,
                 float timeout = -1., bool track_filename = false,
	         size_t buffersize = 1024*1024);
	G3Reader(const std::vector<std::string> &filenames, int n_frames_to_read = -1,
                 float timeout = -1., bool track_filename = false,
	         size_t buffersize = 1024*1024);

	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
	off_t Seek(off_t offset);
	off_t Tell();

private:
	void StartFile(const std::string &path);
	bool prefix_file_;
	std::string cur_file_;
	std::deque<std::string> filename_;
	std::istream stream_;
	int n_frames_to_read_;
	int n_frames_read_;
	int n_frames_cur_;
	float timeout_;
	bool track_filename_;
	size_t buffersize_;

	SET_LOGGER("G3Reader");
};

G3_POINTERS(G3Reader);

#endif
