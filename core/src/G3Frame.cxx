#include <G3Frame.h>
#include <G3Data.h>
#include <G3Quat.h>
#include <serialization.h>
#include <pybindings.h>
#include <container_pybindings.h>

#include <stdlib.h>
#include <cxxabi.h>
#include <algorithm>

extern "C" unsigned long crc32c(unsigned long crc, const uint8_t *buf, unsigned int len);

template <typename T>
static std::string cxx_demangle(const T& v)
{
	int err = 0;
	const char *name = typeid(v).name();
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
	return cxx_demangle(*this);
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
		py::gil_scoped_acquire gil;

		try {
			auto pyobj = py::cast(obj);
			auto cls = pyobj.attr("__class__");

			return cls.attr("__module__").cast<std::string>() +
			    "." + cls.attr("__name__").cast<std::string>();
		} catch (...) {
			// Fall through to C++ name
		}
	}

	return cxx_demangle(*obj);
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
	case G3Frame::Ephemeris:
		ft_str = "Ephemeris";
		break;
	case G3Frame::LightCurve:
		ft_str = "LightCurve";
		break;
	case G3Frame::Statistics:
		ft_str = "Statistics";
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
		// Right-justify the character string in the constant in native
		// endianness.
		for (int j = 24; j >= 0; j -= 8) {
			char c = char((uint32_t(frame_type) >> j) & 0xff);
			if (c == 0)
				continue;
			ft_str += c;
		}
	}

	os << ft_str;
	return os;
}

template <typename T>
void G3Frame::saves(T &os) const
{
	using cereal::make_nvp;

	cereal::PortableBinaryOutputArchive ar(os);
	ar << make_nvp("frame", *this);
}

template <class A>
void G3Frame::save(A &ar, unsigned v) const
{
	using cereal::make_nvp;
	uint32_t crc(0), size(map_.size());

	ar << make_nvp("size", size);
	ar << make_nvp("type", (uint32_t)type);
	for (auto i = map_.begin(); i != map_.end(); i++) {
		// If no encoded frameobject, serialize it
		blob_encode(i->second);

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
void G3Frame::loads(T &is)
{
	using cereal::make_nvp;

	cereal::PortableBinaryInputArchive ar(is);
	ar >> make_nvp("frame", *this);
}

template <class A>
void G3Frame::load(A &ar, unsigned v)
{
	using cereal::make_nvp;
	G3_CHECK_VERSION(v);

	int size;
	uint32_t xtype, testcrc;
	uint32_t crc(0);

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
		blob.blob = std::make_shared<std::vector<char> >();
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

void G3Frame::DropBlobs(bool decode_all) const
{
	for (auto i = map_.begin(); i != map_.end(); i++) {
		if (!i->second.frameobject) {
			// Decode if asked for and possible, otherwise leave
			// the blob in place.
			if (decode_all) {
				try {
					blob_decode(i->second);
				} catch (cereal::Exception &e) {
					continue;
				}
			} else {
				continue;
			}
		}

		i->second.blob.reset();
	}
}

void G3Frame::GenerateBlobs(bool drop_objects) const
{
	for (auto i = map_.begin(); i != map_.end(); i++) {
		// If no encoded frameobject, serialize it
		blob_encode(i->second);

		if (drop_objects)
			i->second.frameobject.reset();
	}
}

void G3Frame::DropObjects() const
{
	for (auto i = map_.begin(); i != map_.end(); i++) {
		if (i->second.blob)
			i->second.frameobject.reset();
	}
}

void G3Frame::blob_decode(struct blob_container &blob)
{
	using cereal::make_nvp;
	G3FrameObjectPtr ptr;

	if (blob.frameobject)
		return;

	G3BufferInputStream is(*blob.blob);
	cereal::PortableBinaryInputArchive item_ar(is);
	item_ar >> make_nvp("val", ptr);
	blob.frameobject = ptr;

	// Drop really big (> 128 MB) blobs at decode time. The savings
	// in CPU is not really worth the huge memory hit.
	if (blob.blob->size() > 128*1024*1024)
		blob.blob.reset();
}

void G3Frame::blob_encode(struct blob_container &blob)
{
	using cereal::make_nvp;

	if (blob.blob)
		return;

	// If no encoded frameobject, serialize it
	blob.blob = std::make_shared<std::vector<char> >();
	G3BufferOutputStream item_os(*blob.blob);
	cereal::PortableBinaryOutputArchive item_ar(item_os);
	item_ar << make_nvp("val", blob.frameobject);
	item_os.flush();
}

template void G3Frame::loads(G3BufferInputStream &);
template void G3Frame::loads(std::istream &);
template void G3Frame::loads(std::istringstream &);

template void G3Frame::saves(G3BufferOutputStream &) const;
template void G3Frame::saves(std::ostream &) const;
template void G3Frame::saves(std::ostringstream &) const;

G3_SPLIT_SERIALIZABLE_CODE(G3Frame);

G3FramePtr
g3frame_char_constructor(std::string max_4_chars)
{
	if (max_4_chars.size() > 4) {
		throw py::value_error("Ad-hoc frame type must be 4 "
		    "or fewer characters.");
	}

	// Right-justify the character string in the constant in native
	// endianness, such that code = max_4_chars[0] for a 1-character
	// code.
	uint32_t code = 0;
	for (int i = max_4_chars.size()-1, j = 0; i >= 0; i--, j += 8)
		code |= (uint32_t(max_4_chars[i]) << j);

	return G3FramePtr(new G3Frame(G3Frame::FrameType(code)));
}

static void g3frame_python_put(G3Frame &f, const std::string &name, const py::object &obj)
{
	if (py::isinstance<py::bool_>(obj))
		f.Put(name, std::make_shared<G3Bool>(obj.cast<bool>()));
	else if (py::isinstance<py::int_>(obj))
		f.Put(name, std::make_shared<G3Int>(obj.cast<int64_t>()));
	else if (py::isinstance<py::float_>(obj))
		f.Put(name, std::make_shared<G3Double>(obj.cast<double>()));
	else if (py::isinstance<Quat>(obj))
		f.Put(name, std::make_shared<G3Quat>(obj.cast<Quat>()));
	else if (py::isinstance<py::str>(obj))
		f.Put(name, std::make_shared<G3String>(obj.cast<std::string>()));
	else
		throw py::type_error(
		    "Object is not a G3FrameObject derivative or a plain-old-data type");
}

static py::object g3frame_python_get(const G3Frame &f, const std::string &name)
{
	// Python doesn't have a concept of const. Add subterfuge.
	G3FrameObjectConstPtr element = f[name];
	if (!element)
		throw py::key_error(name);

	if (!!std::dynamic_pointer_cast<const G3Int>(element))
		return py::int_(std::dynamic_pointer_cast<const G3Int>(element)->value);
	else if (!!std::dynamic_pointer_cast<const G3Double>(element))
		return py::float_(std::dynamic_pointer_cast<const G3Double>(element)->value);
	else if (!!std::dynamic_pointer_cast<const G3String>(element))
		return py::str(std::dynamic_pointer_cast<const G3String>(element)->value);
	else if (!!std::dynamic_pointer_cast<const G3Bool>(element))
		return py::bool_(std::dynamic_pointer_cast<const G3Bool>(element)->value);
	else if (!!std::dynamic_pointer_cast<const G3Quat>(element))
		return py::cast(std::dynamic_pointer_cast<const G3Quat>(element)->value);
	else
		return py::cast(element);
}

static py::object
g3frame_python_get_default(const G3Frame &f, const std::string &name, const py::object &d)
{
	try {
		return g3frame_python_get(f, name);
	} catch (py::key_error &) {
		return d;
	}
}

static py::object
g3frame_python_pop(G3Frame &f, const std::string &name)
{
	auto v = g3frame_python_get(f, name);
	f.Delete(name);
	return v;
}

static py::object
g3frame_python_pop_default(G3Frame &f, const std::string &name, const py::object &d)
{
	auto v = g3frame_python_get_default(f, name, d);
	f.Delete(name);
	return v;
}

static std::string g3frame_str(const G3Frame &f)
{
	std::ostringstream oss;
	oss << f;
	return oss.str();
}

static py::bytes g3frame_hash(const py::object &obj)
{
	return obj.attr("__getstate__")().cast<py::tuple>()[1];
}


// Iterators for keys/views/items
class PyFrameIterator {
public:
	using value_type = std::pair<const std::string, py::object>;

	PyFrameIterator(const G3Frame &frame, const G3Frame::key_iterator &it) :
	    frame(frame), it(it) {}

	bool operator ==(const PyFrameIterator & other) const { return it == other.it; }
	PyFrameIterator operator++() { ++it; return *this; }
	value_type operator*() const { return {*it, g3frame_python_get(frame, *it)}; }

private:
	const G3Frame &frame;
	G3Frame::key_iterator it;
};

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
PYBIND11_NAMESPACE_BEGIN(detail)

template <>
struct KeysViewImpl<G3Frame> : public py::detail::keys_view {
	explicit KeysViewImpl(G3Frame &frame) : frame(frame) {}
	size_t len() override { return frame.size(); }
	iterator iter() override { return make_iterator(frame.begin(), frame.end()); }
	bool contains(const py::handle &k) override {
		try {
			return frame.Has(k.cast<std::string>());
		} catch (const py::cast_error &) {
			return false;
		}
	}
	G3Frame &frame;
};

template <>
struct ValuesViewImpl<G3Frame> : public py::detail::values_view {
	explicit ValuesViewImpl(G3Frame &frame) : frame(frame) {}
	size_t len() override { return frame.size(); }
	iterator iter() override {
		return make_value_iterator(PyFrameIterator(frame, frame.begin()),
		    PyFrameIterator(frame, frame.end()));
	}
	G3Frame &frame;
};

template <>
struct ItemsViewImpl<G3Frame> : public py::detail::items_view {
	explicit ItemsViewImpl(G3Frame &frame) : frame(frame) {}
	size_t len() override { return frame.size(); }
	iterator iter() override {
		return make_iterator(PyFrameIterator(frame, frame.begin()),
		    PyFrameIterator(frame, frame.end()));
	}
	G3Frame &frame;
};

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)


PYBINDINGS("core", scope) {
	// Internal stuff
	register_class<G3FrameObject>(scope, "G3FrameObject",
	  "Base class for objects that can be added to a frame. All such "
	  "must inherit from G3FrameObject in C++. Pickle hooks are overridden "
	  "to use the fast internal serialization")
	    .def("Description", &G3FrameObject::Description,
	      "Long-form human-readable description of the object")
	    .def("Summary", &G3FrameObject::Summary,
	      "Short (one-line) description of the object")
	    .def("__str__", &G3FrameObject::Summary)
	    .def(g3frameobject_picklesuite<G3FrameObject>())
	    .def_property_readonly("hash", &g3frame_hash,
	      "Return the serialized representation of the object")
	;

	register_enum<G3Frame::FrameType, G3Frame::None>(scope, "G3FrameType")
	    .value("Timepoint",       G3Frame::Timepoint)
	    .value("Housekeeping",    G3Frame::Housekeeping)
	    .value("Observation",     G3Frame::Observation)
	    .value("Scan",            G3Frame::Scan)
	    .value("Map",             G3Frame::Map)
	    .value("InstrumentStatus",G3Frame::InstrumentStatus)
	    .value("PipelineInfo",    G3Frame::PipelineInfo)
	    .value("EndProcessing",   G3Frame::EndProcessing)
	    .value("Calibration",     G3Frame::Calibration)
	    .value("Wiring",          G3Frame::Wiring)
	    .value("GcpSlow",         G3Frame::GcpSlow)
	    .value("Ephemeris",       G3Frame::Ephemeris)
	    .value("LightCurve",      G3Frame::LightCurve)
	    .value("Statistics",      G3Frame::Statistics)
	    .value("none",            G3Frame::None)
	    .def_property_readonly("key",
	        [](G3Frame::FrameType &v) { return static_cast<char>(v); },
	        "Return G3FrameType's C++ enum character key")
	;
	register_vector_of<G3Frame::FrameType>(scope, "FrameType");

	auto cls = register_class<G3Frame>(scope, "G3Frame",
	  "Frames are the core datatype of the analysis software. They behave "
	  "like Python dictionaries except that they can only store subclasses "
	  "of core.G3FrameObject and the dictionary keys must be strings. "
	  "Pickling and unpickling them uses internal serialization and is "
	  "extremely fast."
	  "\n\n"
	  "In addition to dictionary-like contents, frames have a type code "
	  "(G3Frame.Type) that designates what kind of data are contained in "
	  "it. These types usually indicates information that changes at "
	  "different rates.")
	    .def(py::init<>())
	    .def(py::init<G3Frame::FrameType>())
	    .def(py::init<G3Frame>())
	    .def(py::init(&g3frame_char_constructor), py::arg("adhoctypecode"),
	      "Create a frame with an ad-hoc (non-standard) type code. "
	      "Use sparingly and with care.")
	    .def_readwrite("type", &G3Frame::type,
	      "Type code for frame. See general G3Frame docstring.")
	    .def_readonly("_filename", &G3Frame::_filename,
	      "Source filename for frame, if read in using G3Reader. This attribute is "
	      "fragile, use at your own risk.")
	    .def("__setitem__", &G3Frame::Put)
	    .def("__setitem__", &g3frame_python_put)
	    .def("__getitem__", &g3frame_python_get)
	    .def("get", &g3frame_python_get_default, py::arg("key"),
	      py::arg("default")=py::none())
	    .def("__delitem__", &G3Frame::Delete)
	    .def("pop", &g3frame_python_pop)
	    .def("pop", &g3frame_python_pop_default)
	    .def("__contains__", (bool (G3Frame::*)(const std::string &) const)
	       &G3Frame::Has)
	    .def("__str__", &g3frame_str)
	    .def("__len__", &G3Frame::size)
	    .def("__iter__",
	       [](G3Frame &f) { return py::make_iterator(f.begin(), f.end()); },
	       py::keep_alive<0, 1>())
	    .def("drop_blobs", &G3Frame::DropBlobs, py::arg("decode_all")=false,
	      "Drop all serialized data either for already-decoded objects (default) "
	      "or all objects after decoding them (if decode_all is true). "
	      "Saves memory at the expense of CPU time if reserialized.")
	    .def("generate_blobs", &G3Frame::GenerateBlobs, py::arg("drop_objects")=false,
	      "Force immediate serialization of all objects. Will save some "
	      "CPU time later during serialization of the frame in exchange "
	      "for spending the exact same amount of CPU time right now.")
	    .def("drop_objects", &G3Frame::DropObjects,
	      "Drop all decoded objects in favor of their serialized copies, "
	      "where those serialized copies already exist. Saves memory for "
	      "frames about to be written at the expense of CPU time to "
	      "re-decode them if they are accessed again later.")
	    .def(g3frameobject_picklesuite<G3Frame>())
	    .def_property_readonly("hash", &g3frame_hash,
	      "Return the serialized representation of the frame")
	;
	register_map_iterators<G3Frame>(scope, cls);
	register_vector_of<G3FramePtr>(scope, "Frame");
	register_vector_of<G3FrameObjectPtr>(scope, "FrameObject");
}
