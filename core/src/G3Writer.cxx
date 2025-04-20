#include <pybindings.h>
#include <G3Writer.h>
#include <dataio.h>

G3Writer::G3Writer(std::string filename,
    std::vector<G3Frame::FrameType> streams,
    bool append, size_t buffersize) :
    filename_(filename), stream_(nullptr), streams_(streams)
{
	g3_ostream_to_path(stream_, filename, append, buffersize);
}

void G3Writer::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
{
	// If in Python context, release the GIL while writing out.
	// Note we must force serialization of the frames, first,
	// because in some cases serialization requires calling back
	// out to Python.
	frame->GenerateBlobs();

	G3PythonContext ctx("G3Writer", false);

	if (frame->type == G3Frame::EndProcessing)
		stream_.flush();
	else if (streams_.size() == 0 ||
	    std::find(streams_.begin(), streams_.end(), frame->type) !=
	    streams_.end())
		frame->saves(stream_);

	out.push_back(frame);
}

void G3Writer::Flush()
{
	stream_.flush();
}

off_t G3Writer::Tell() {
	return stream_.tellp();
}

PYBINDINGS("core") {
	// Instead of EXPORT_G3MODULE since there is an extra Flush function
	py::class_<G3Writer, py::bases<G3Module>, std::shared_ptr<G3Writer>,
	    boost::noncopyable>("G3Writer",
	      "Writes frames to disk. Frames will be written to the file specified by "
	      "filename. If filename ends in .gz, output will be compressed using gzip. "
	      "To write only some types of frames, pass a list of the desired frame "
	      "types to the second optional argument (streams). If no streams argument "
	      "is given, writes all types of frames. If append is set to True, will "
	      "append frames to its output file rather than overwriting it.",
	py::init<std::string, std::vector<G3Frame::FrameType>, bool, size_t>((py::arg("filename"),
	    py::arg("streams")=std::vector<G3Frame::FrameType>(), py::arg("append")=false,
	    py::arg("buffersize")=1024*1024)))
	.def("flush", &G3Writer::Flush)
	.def("tell", &G3Writer::Tell)
	.def_readonly("__g3module__", true)
	;
}
