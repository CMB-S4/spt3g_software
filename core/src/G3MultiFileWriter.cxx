#include <pybindings.h>
#include <string>
#include <G3Module.h>

#include <dataio.h>

class G3MultiFileWriter : public G3Module {
public:
	G3MultiFileWriter(boost::python::object filename,
	    size_t size_limit,
	    boost::python::object divide_on = boost::python::object());
	~G3MultiFileWriter();
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
	std::string CurrentFile() { return current_filename_; }
private:
	bool CheckNewFile(G3FramePtr frame);

	std::string filename_;
	boost::python::object filename_callback_;
	std::string current_filename_;
	size_t size_limit_;

	std::vector<G3Frame::FrameType> always_break_on_;
	boost::python::object newfile_callback_;

	std::ostream stream_;
	std::vector<G3FramePtr> metadata_cache_;
	int seqno;

	SET_LOGGER("G3MultiFileWriter");
};

G3MultiFileWriter::G3MultiFileWriter(boost::python::object filename,
    size_t size_limit, boost::python::object divide_on)
    : size_limit_(size_limit), stream_(nullptr), seqno(0)
{
	boost::python::extract<std::string> fstr(filename);

	if (fstr.check()) {
		filename_ = fstr();
		g3_check_output_path(filename_);

		if (snprintf(NULL, 0, filename_.c_str(), 0) < 0)
			log_fatal("Cannot format filename. Should be "
			    "outfile-%%03u.g3");
	} else if (PyCallable_Check(filename.ptr())) {
		filename_ = "";
		filename_callback_ = filename;
	} else {
		log_fatal("filename must be either a string with a format "
		    "character for file number or a Python callable that "
		    "returns a string with the signature f(frame, seqno)");
	}

	if (size_limit_ == 0)
		log_fatal("File size limit must be greater than zero");

	if (divide_on.ptr() != Py_None) {
		boost::python::extract<std::vector<G3Frame::FrameType> >
		    type_list_ext(divide_on);

		if (type_list_ext.check())
			always_break_on_ = type_list_ext();
		else if (PyCallable_Check(divide_on.ptr()))
			newfile_callback_ = divide_on;
		else
			log_fatal("divide_on must be either an iterable of "
			    "frame types on which to start a new file "
			    "(e.g. [core.G3FrameType.Observation]) or "
			    "a callable that inspects a frame and returns "
			    "True if a new file should be started and False "
			    "otherwise.");
	}
}

G3MultiFileWriter::~G3MultiFileWriter()
{
	g3_ostream_close(stream_);
}

bool
G3MultiFileWriter::CheckNewFile(G3FramePtr frame)
{
	// If we are already saving data, check file size. Otherwise, open
	// a new file unconditionally.
	if (stream_) {
		bool start_new_ = false;

		if ((size_t)stream_.tellp() > size_limit_)
			start_new_ = true;

		if (newfile_callback_.ptr() != Py_None &&
		    boost::python::extract<bool>(newfile_callback_(frame))())
			start_new_ = true;

		if (std::find(always_break_on_.begin(), always_break_on_.end(),
		    frame->type) != always_break_on_.end())
			start_new_ = true;

		if (!start_new_)
			return false;
	}

	g3_ostream_close(stream_);

	std::string filename;
	if (filename_ != "") {
		int sz = snprintf(NULL, 0, filename_.c_str(), seqno);
		if (sz < 0)
			log_fatal("Cannot format filename. Should be "
			    "outfile-%%03u.g3");
		char *msg = new char[sz + 1];
		snprintf(msg, sz + 1, filename_.c_str(), seqno);
		filename = std::string(msg);
		delete [] msg;
		seqno++;
	} else {
		filename = boost::python::extract<std::string>(
		    filename_callback_(frame, seqno++))();

		g3_check_output_path(filename);
	}

	current_filename_ = filename;
	g3_ostream_to_path(stream_, filename, false);

	for (auto i = metadata_cache_.begin(); i != metadata_cache_.end(); i++)
		(*i)->saves(stream_);

	return true;
}

void G3MultiFileWriter::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
{
	bool new_file(false), meta_cached(false);

	if (frame->type == G3Frame::EndProcessing) {
		g3_ostream_close(stream_);
		goto done;
	}

	// Update metadata cache so new files have the right headers
	if (frame->type != G3Frame::Timepoint && frame->type != G3Frame::Scan) {
		auto i = metadata_cache_.begin();
		while (i != metadata_cache_.end()) {
			if ((*i)->type == frame->type) {
				*i = frame;
				break;
			}
			i++;
		}
		if (i == metadata_cache_.end()) // Nothing yet of this type
			metadata_cache_.push_back(frame);

		meta_cached = true;
	}

	// See if we need a new file (this implicitly opens one if so)
	new_file = CheckNewFile(frame);

	// And out to disk, making sure not to write again if it just went into
	// the metadata cache and onto disk in CheckNewFile()
	if (!new_file || !meta_cached)
		frame->saves(stream_);

done:
	out.push_back(frame);
}

PYBINDINGS("core") {
	using namespace boost::python;

	class_<G3MultiFileWriter, bases<G3Module>, std::shared_ptr<G3MultiFileWriter>,
	    boost::noncopyable>("G3MultiFileWriter",
	      "Writes frames to disk into a sequence of files. Once a file exceeds "
	      "the number of bytes specified in size_limit, it will start a new file. "
	      "Files are named based on filename. If passed a string for filename "
	      "with a printf-style specifier, that specifier will be replaced by a "
	      "zero-indexed sequence number. For example, outfile-%03u.g3.gz would "
	      "produce a sequence of files named outfile-000.g3.gz, outfile-001.g3.gz, "
	      "etc. Alternatively, you can pass a callable that is passed the first "
	      "frame in the new file and the sequence number and returns a path to "
	      "the new file. Any frames besides Timepoint and Scan frames have the "
	      "most recent frame of each type prepended to all new files.\n\n"
	      "More complex behavior can be obtained with the optional divide_on "
	      "argument. This can be an iterable of frame types (e.g. "
	      "[core.G3FrameType.Observation]) or a callable. In the iterable case, "
	      "the presence of any frame with a type in the list will cause the "
	      "creation of a new file even if the file size threshold has not yet "
	      "been met. This is useful to create files based on, for example, "
	      "observation boundaries. For more flexibility, you can also pass a "
	      "python callable as divide_on. This callable will be passed each "
	      "frame in turn. If it returns True (or something with positive "
	      "truth-value), a new file will be started at that frame.",
	init<object, size_t, optional<object> >((arg("filename"),
	    arg("size_limit"), arg("divide_on")=object())))
	.def_readonly("current_file", &G3MultiFileWriter::CurrentFile)
	.def_readonly("__g3module__", true)
	;
}
