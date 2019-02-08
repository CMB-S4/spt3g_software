#include <G3Frame.h>
#include <serialization.h>
#include <pybindings.h>

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <sstream>
#include <stdlib.h>
#include <cxxabi.h>
#include <algorithm>

extern "C" unsigned long crc32c(unsigned long crc, const uint8_t *buf, unsigned int len);

// Work around crummy slow implementation of xsgetn() in libc++
class fast_streambuf : public std::basic_streambuf<char>
{
public:
	fast_streambuf(std::vector<char> &vec) : basic_streambuf() {
		setg(&vec[0], &vec[0], &vec[vec.size()]);
	}
protected:
	std::streamsize xsgetn(char_type *buf, std::streamsize n) {
		std::streamsize to_read = egptr() - gptr();
		if (to_read > n)
			to_read = n;
		memcpy(buf, gptr(), to_read);
		gbump(to_read);
		return to_read;
	}
};

static std::string cxx_demangle(const char *name)
{
	int err = 0;
	char *demangled;

	demangled = abi::__cxa_demangle(name, NULL, NULL, &err);

	std::string demangled_name((err == 0) ? demangled : name);
	free(demangled);

	return demangled_name;
}

std::ostream& operator<<(std::ostream& os, const G3FrameObject &frame)
{
	return (os << frame.Description());
}

std::string G3FrameObject::Description() const
{
	return cxx_demangle(typeid(*this).name());
}

std::string G3FrameObject::Summary() const
{
	return Description();
}

template <class A> void G3FrameObject::serialize(A &ar, const unsigned v)
{
}

G3_SERIALIZABLE_CODE(G3FrameObject);

G3Frame::G3Frame(G3Frame::FrameType t) : type(t)
{}

G3FrameObjectConstPtr G3Frame::operator [](const std::string &name) const
{
	G3MapType::iterator iter;

	iter = map_.find(name);

	if (iter == map_.end())
		return G3FrameObjectConstPtr();

	blob_decode(iter->second);
	return iter->second.frameobject;
}

void G3Frame::Put(const std::string &name, G3FrameObjectConstPtr value)
{
	std::pair<G3MapType::iterator, bool> result;
	struct blob_container blob;
	blob.frameobject = value;

	if (!value)
		log_fatal("Cannot add None to frame");

	result = map_.insert(G3MapType::value_type(name, blob));

	if (!result.second)
		log_fatal("Previously existing key \"%s\"", name.c_str());
}

void G3Frame::Delete(const std::string &name)
{

	map_.erase(name);
}

bool G3Frame::Has(const std::string &name) const
{
	return map_.find(name) != map_.end();
}

std::vector<std::string> G3Frame::Keys() const
{
	std::vector<std::string> keys;

	for (auto i = map_.begin(); i != map_.end(); i++)
		keys.push_back(i->first);

	return keys;
}

G3Frame &G3Frame::operator = (const G3Frame &copy)
{
	map_ = copy.map_;
	type = copy.type;
	return *this;
}

static std::string FrameObjectClassName(G3FrameObjectConstPtr obj)
{
	// Try to give python name if possible
	if (Py_IsInitialized()) {
		try {
			boost::python::object pyobj(
			    boost::const_pointer_cast<G3FrameObject>(obj));

			return
			    boost::python::extract<std::string>(
			     pyobj.attr("__class__").attr("__module__"))() +
			     "." +
			     boost::python::extract<std::string>(
			      pyobj.attr("__class__").attr("__name__"))();
		} catch (...) {
			// Fall through to C++ name
		}
	}

	return cxx_demangle(typeid(*obj).name());
}

std::ostream& operator<<(std::ostream& os, const G3Frame &frame)
{
	// sorts the keys for printing
	std::vector<std::string> key_vec(frame.map_.size());
	size_t j = 0;
	for (auto it = frame.map_.begin(); it != frame.map_.end(); it++) {
		key_vec[j] = it->first;
		j++;
	}
	std::sort(key_vec.begin(), key_vec.end());

	os << "Frame (" << frame.type << ") [" << std::endl;
	for (auto i = key_vec.begin(); i != key_vec.end(); i++) {
		try {
			G3Frame::blob_decode(frame.map_.at(*i));
			G3FrameObjectConstPtr obj = 
			    frame.map_.at(*i).frameobject;
			os << "\"" << *i << "\" (" <<
			    FrameObjectClassName(obj) <<
			    ") => " << obj->Summary() << std::endl;
		} catch (cereal::Exception &e) {
			// Cereal provides really verbose what(). Truncate
			// to only the first line.
			std::string what(e.what());

			os << "\"" << *i << "\" (unreadable, " <<
			    frame.map_.at(*i).blob->size() << " bytes) => " <<
			    what.substr(0, what.find('\n')) << std::endl;
		}
	}
	os << "]";
	return os;
}

std::ostream& operator<<(std::ostream& os, const G3Frame::FrameType &frame_type)
{
	std::string ft_str;

	switch (frame_type) {
	case G3Frame::Timepoint:
		ft_str = "Timepoint";
		break;
	case G3Frame::Housekeeping:
		ft_str = "Housekeeping";
		break;
	case G3Frame::Observation:
		ft_str = "Observation";
		break;
	case G3Frame::Scan:
		ft_str = "Scan";
		break;
	case G3Frame::Map:
		ft_str = "Map";
		break;
	case G3Frame::InstrumentStatus:
		ft_str = "InstrumentStatus";
		break;
	case G3Frame::PipelineInfo:
		ft_str = "PipelineInfo";
		break;
	case G3Frame::GcpSlow:
		ft_str = "GcpSlow";
		break;
	case G3Frame::EndProcessing:
	        ft_str = "EndProcessing";
		break;
	case G3Frame::Calibration:
		ft_str = "Calibration";
		break;
	case G3Frame::Wiring:
		ft_str = "Wiring";
		break;
	case G3Frame::None:
		ft_str = "None";
		break;
	default:
		ft_str = frame_type;
	}

	os << ft_str;
	return os;
}

template <typename T>
void G3Frame::save(T &os) const
{
	using cereal::make_nvp;
	uint32_t crc = 0;
	uint32_t version(1), size(map_.size());

	cereal::PortableBinaryOutputArchive ar(os);
	ar << make_nvp("version", version);
	ar << make_nvp("size", size);
	ar << make_nvp("type", (uint32_t)type);
	for (auto i = map_.begin(); i != map_.end(); i++) {
		if (!i->second.blob) {
			// If no encoded frameobject, serialize it
			i->second.blob = boost::make_shared<std::vector<char> >();
			boost::iostreams::stream<
			  boost::iostreams::back_insert_device<std::vector<char> > >
			  item_os(*i->second.blob);
			cereal::PortableBinaryOutputArchive item_ar(item_os);
			item_ar << make_nvp("val", i->second.frameobject);
			item_os.flush();
		}
		ar << make_nvp("name", i->first);
		crc = crc32c(crc, (const uint8_t *)i->first.c_str(),
		    i->first.size());
		ar << make_nvp("blob", *i->second.blob);
		crc = crc32c(crc, (const uint8_t *)&(*i->second.blob)[0],
		    i->second.blob->size());
	}
	ar << make_nvp("crc", crc);
}

template <typename T>
void G3Frame::load(T &is)
{
	using cereal::make_nvp;

	cereal::PortableBinaryInputArchive ar(is);
	int version, size;
	uint32_t xtype;
	uint32_t crc(0), testcrc;

	ar >> make_nvp("version", version);
	ar >> make_nvp("size", size);
	ar >> make_nvp("type", xtype);
	type = FrameType(xtype);

	map_.clear();
	for (int i = 0; i < size; i++) {
		struct blob_container blob;
		std::string name;
		G3FrameObjectPtr ptr;

		ar >> make_nvp("name", name);
		crc = crc32c(crc, (const uint8_t *)name.c_str(), name.size());
		blob.blob = boost::make_shared<std::vector<char> >();
		ar >> make_nvp("blob", *blob.blob);
		crc = crc32c(crc, (const uint8_t *)&(*blob.blob)[0],
		    blob.blob->size());
		map_.insert(G3MapType::value_type(name, blob));
	}
	ar >> make_nvp("crc", testcrc);

	if (testcrc != crc)
		log_fatal("Recorded CRC (%#x) does not match calculated (%#x)",
		    testcrc, crc);
}

template <>
void G3Frame::load(std::vector<char> &data)
{
	fast_streambuf sb(data);
	std::istream is(&sb);

	load(is);
}

void G3Frame::blob_decode(struct blob_container &blob)
{
	using cereal::make_nvp;
	G3FrameObjectPtr ptr;

	if (blob.frameobject)
		return;

	fast_streambuf sb(*blob.blob);
	std::istream is(&sb);
	cereal::PortableBinaryInputArchive item_ar(is);
	item_ar >> make_nvp("val", ptr);
	blob.frameobject = ptr;

	// Drop really big (> 128 MB) blobs at decode time. The savings
	// in CPU is not really worth the huge memory hit.
	if (blob.blob->size() > 128*1024*1024)
		blob.blob.reset();
}

template void G3Frame::load(boost::iostreams::filtering_istream &);
template void G3Frame::load(std::istream &);
template void G3Frame::load(std::istringstream &);

template void G3Frame::save(boost::iostreams::filtering_ostream &) const;
template void G3Frame::save(std::ostream &) const;
template void G3Frame::save(std::ostringstream &) const;

