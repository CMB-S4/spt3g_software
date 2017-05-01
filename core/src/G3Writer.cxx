#include <pybindings.h>
#include <G3Writer.h>

#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/filesystem.hpp>

G3Writer::G3Writer(std::string filename) : filename_(filename)
{
	boost::filesystem::path fpath(filename);
	if ( fpath.empty() || ( fpath.has_parent_path() && ! boost::filesystem::exists( fpath.parent_path() )))
		throw std::runtime_error(std::string("Parent path does not exist: ") + fpath.parent_path().string());

	if (boost::algorithm::ends_with(filename, ".gz"))
		stream_.push(boost::iostreams::gzip_compressor());
	stream_.push(boost::iostreams::file_sink(filename, std::ios::binary));
}

G3Writer::~G3Writer()
{
	stream_.reset();
}

void G3Writer::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
{
	if (frame->type == G3Frame::EndProcessing)
		stream_.reset();
	else
		frame->save(stream_);
	out.push_back(frame);
}

EXPORT_G3MODULE("core", G3Writer, init<std::string>(args("filename")),
    "Writes frames to disk");

