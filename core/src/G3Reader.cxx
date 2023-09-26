#include <pybindings.h>
#include <dataio.h>
#include <G3Reader.h>

#include <boost/filesystem.hpp>

G3Reader::G3Reader(std::string filename, int n_frames_to_read,
    float timeout, bool track_filename) :
    prefix_file_(false), n_frames_to_read_(n_frames_to_read),
    n_frames_read_(0), timeout_(timeout), track_filename_(track_filename)
{
	boost::filesystem::path fpath(filename);
	if (filename.find("://") == std::string::npos &&
	   (!boost::filesystem::exists(fpath) ||
	    !boost::filesystem::is_regular_file(fpath)))
		log_fatal("Could not find file %s", filename.c_str());
	StartFile(filename);
}

G3Reader::G3Reader(std::vector<std::string> filename, int n_frames_to_read,
    float timeout, bool track_filename) :
    prefix_file_(false), n_frames_to_read_(n_frames_to_read),
    n_frames_read_(0), timeout_(timeout), track_filename_(track_filename)
{
	if (filename.size() == 0)
		log_fatal("Empty file list provided to G3Reader");

	for (auto i = filename.begin(); i != filename.end(); i++){
		boost::filesystem::path fpath(*i);
		if (i->find("://") == std::string::npos &&
		   (!boost::filesystem::exists(fpath) ||
		    !boost::filesystem::is_regular_file(fpath)))
			log_fatal("Could not find file %s", i->c_str());
		filename_.push_back(*i);
	}
	StartFile(filename_.front());
	filename_.pop_front();
}

void G3Reader::StartFile(std::string path)
{
	log_info("Starting file %s\n", path.c_str());
	cur_file_ = path;
	(void) g3_istream_from_path(stream_, path, timeout_);
}

void G3Reader::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
{
	// Bail (ends processing) if too many frames have passed by
	if (!frame && n_frames_to_read_ > 0 && 
	    n_frames_read_ >= n_frames_to_read_)
		return;

	// If this is not a driving module and we haven't already, prepend the
	// file to the input frame
	if (frame) {
		if (!prefix_file_) {
			prefix_file_ = true;
			std::deque<G3FramePtr> subdeque;
			while (true) {
				Process(G3FramePtr(), subdeque);
				if (subdeque.size() == 0)
					break;
				for (auto i = subdeque.begin();
				    i != subdeque.end(); i++)
					out.push_back(*i);
				subdeque.clear();
			}
		}
		out.push_back(frame);
	}

	// For python threaded contexts, release the GIL while
	// blocking for data (which can occur on stream.peek as well
	// as load(stream_)).  Only do this if we're actually in a
	// Python interpreter (pure C++ applications will fail the
	// Py_IsInitialized test).  Make sure all paths out of this
	// function reacquire the lock, if it was released.
	PyThreadState *_save = nullptr;
	if (Py_IsInitialized())
		_save = PyEval_SaveThread();

	if (stream_.peek() == EOF) {
		if (filename_.size() > 0) {
			StartFile(filename_.front());
			filename_.pop_front();
		} else {
			// Stop processing
			if (_save != nullptr)
				PyEval_RestoreThread(_save);
			return;
		}
	}
	frame = G3FramePtr(new G3Frame);
	try {
		frame->load(stream_);
	} catch (...) {
		log_error("Exception raised while reading file %s",
		    cur_file_.c_str());
		if (_save != nullptr)
			PyEval_RestoreThread(_save);
		throw;
	}
	if (_save != nullptr)
		PyEval_RestoreThread(_save);

	if (track_filename_)
		frame->_filename = cur_file_;
	out.push_back(frame);

	n_frames_read_++;
}

off_t G3Reader::Seek(off_t offset) {
	return boost::iostreams::seek(stream_, offset, std::ios_base::beg);
}

off_t G3Reader::Tell() {
	return boost::iostreams::seek(stream_, 0, std::ios_base::cur);
}

PYBINDINGS("core") {
	using namespace boost::python;

	// Instead of EXPORT_G3MODULE since there are two constructors
	class_<G3Reader, bases<G3Module>, boost::shared_ptr<G3Reader>,
	    boost::noncopyable>("G3Reader",
	      "Read frames from disk. Takes either the path to a file to read "
	      "or an iterable of files to be read in sequence. If "
	      "n_frames_to_read is greater than zero, will stop after "
	      "n_frames_to_read frames rather than at the end of the file[s]. "
	      "The timeout parameter can used to enable socket timeout for tcp "
	      "streams, resulting in EOF behavior on expiry; unfortunately this "
	      "cannot be used for polling, you have to close the connection. "
	      "Use the `tell` and `seek` methods to record the position of and "
	      "seek to the beginning of a particular frame in the file.  Set "
	      "track_filename to True to record the filename for each frame in "
	      "the ._filename attribute (fragile).",
	init<std::string, int, float, bool>((arg("filename"),
	    arg("n_frames_to_read")=0,arg("timeout")=-1.,arg("track_filename")=false)))
	.def(init<std::vector<std::string>, int, float, bool>((arg("filename"),
	    arg("n_frames_to_read")=0, arg("timeout")=-1.,arg("track_filename")=false)))
	.def("tell", &G3Reader::Tell)
	.def("seek", &G3Reader::Seek)
	.def_readonly("__g3module__", true)
	;
}

