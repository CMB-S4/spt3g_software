#include <pybindings.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <string>
#include <G3Module.h>

#include "counter64.hpp"

class G3MultiFileWriter : public G3Module {
public:
	G3MultiFileWriter(boost::python::object filename,
	    size_t size_limit,
	    boost::python::object divide_on = boost::python::object());
	virtual ~G3MultiFileWriter();
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
private:
	bool CheckNewFile(G3FramePtr frame);

	std::string filename_;
	boost::python::object filename_callback_;
	size_t size_limit_;

	std::vector<G3Frame::FrameType> always_break_on_;
	boost::python::object newfile_callback_;

	boost::iostreams::filtering_ostream stream_;
	std::vector<G3FramePtr> metadata_cache_;
	int seqno;

	SET_LOGGER("G3MultiFileWriter");
};

G3MultiFileWriter::G3MultiFileWriter(boost::python::object filename,
    size_t size_limit, boost::python::object divide_on)
    : size_limit_(size_limit), seqno(0)
{
	boost::python::extract<std::string> fstr(filename);

	if (fstr.check()) {
		filename_ = fstr();
		boost::filesystem::path fpath(filename_);
		if (fpath.empty() || (fpath.has_parent_path() &&
		    !boost::filesystem::exists(fpath.parent_path())))
			log_fatal("Parent path does not exist: %s",
			    fpath.parent_path().string().c_str());

		boost::format f(filename_);
		try {
			f % 0;
		} catch (const std::exception &e) {
			log_fatal("Cannot format filename. Should be "
			    "outfile-%%03u.g3");
		}
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

bool
G3MultiFileWriter::CheckNewFile(G3FramePtr frame)
{
	// If we are already saving data, check file size. Otherwise, open
	// a new file unconditionally.
	if (!stream_.empty()) {
		bool start_new_ = false;

		boost::iostreams::counter64 *counter =
		    stream_.component<boost::iostreams::counter64>(
		    stream_.size() - 2);
		if (!counter)
			log_fatal("Could not get stream counter");

		if (counter->characters() > size_limit_)
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

	stream_.reset();

	std::string filename;
	if (filename_ != "") {
		boost::format f(filename_);
		try {
			f % seqno++;
		} catch (const std::exception &e) {
			log_fatal("Cannot format filename. Should be "
			    "outfile-%%03u.g3");
		}
		filename = f.str();
	} else {
		filename = boost::python::extract<std::string>(
		    filename_callback_(frame, seqno++))();

		boost::filesystem::path fpath(filename);
		if (fpath.empty() || (fpath.has_parent_path() &&
		    !boost::filesystem::exists(fpath.parent_path())))
			log_fatal("Parent path does not exist: %s",
			    fpath.parent_path().string().c_str());

	}

	if (boost::algorithm::ends_with(filename, ".gz"))
		stream_.push(boost::iostreams::gzip_compressor());
	stream_.push(boost::iostreams::counter64());
	stream_.push(boost::iostreams::file_sink(filename, std::ios::binary));

	for (auto i = metadata_cache_.begin(); i != metadata_cache_.end(); i++)
		(*i)->save(stream_);

	return true;
}

G3MultiFileWriter::~G3MultiFileWriter()
{
	stream_.reset();
}

void G3MultiFileWriter::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
{
	bool new_file(false), meta_cached(false);

	if (frame->type == G3Frame::EndProcessing) {
		stream_.reset();
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
		frame->save(stream_);

done:
	out.push_back(frame);
}

EXPORT_G3MODULE("core", G3MultiFileWriter, (init<boost::python::object, size_t, optional<boost::python::object> >(args("filename", "size_limit", "divide_on"))),
    "Writes frames to disk into a sequence of files. Once a file exceeds the number of bytes specified in size_limit, it will start a new file. Files are named based on filename. If passed a string for filename with a printf-style specifier, that specifier will be replaced by a zero-indexed sequence number. For example, outfile-%03u.g3.gz would produce a sequence of files named outfile-000.g3.gz, outfile-001.g3.gz, etc. Alternatively, you can pass a callable that is passed the first frame in the new file and the sequence number and returns a path to the new file. Any frames besides Timepoint and Scan frames have the most recent frame of each type prepended to all new files.\n\n"
    "More complex behavior can be obtained with the optional divide_on argument. This can be an iterable of frame types (e.g. [core.G3FrameType.Observation]) or a callable. In the iterable case, the presence of any frame with a type in the list will cause the creation of a new file even if the file size threshold has not yet been met. This is useful to create files based on, for example, observation boundaries. For more flexibility, you can also pass a python callable as divide_on. This callable will be passed each frame in turn. If it returns True (or something with positive truth-value), a new file will be started at that frame."
);

