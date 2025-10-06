#ifndef _G3_SERIALIZATION_H
#define _G3_SERIALIZATION_H

#include <sstream>

#include <cereal/archives/portable_binary.hpp>

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/complex.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>

#include <pybindings.h>
#include <G3Logging.h>

#define G3_SERIALIZABLE_CODE(x) \
template void x::serialize(cereal::PortableBinaryOutputArchive &, unsigned); \
template void x::serialize(cereal::PortableBinaryInputArchive &, unsigned); \

#define G3_SPLIT_SERIALIZABLE_CODE(x) \
template void x::save(cereal::PortableBinaryOutputArchive &, unsigned) const; \
template void x::load(cereal::PortableBinaryInputArchive &, unsigned); \

// Work around crummy slow implementation of xsgetn() in libc++
class G3InputStreamBuffer : public std::basic_streambuf<char>
{
public:
	G3InputStreamBuffer(std::vector<char> &vec) : basic_streambuf() {
		setg(vec.data(), vec.data(), vec.data() + vec.size());
	}
	G3InputStreamBuffer(char *buf, size_t len) : basic_streambuf() {
		setg(buf, buf, buf + len);
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

class G3BufferInputStream : public std::istream
{
public:
	G3BufferInputStream(std::vector<char> &vec) : std::istream(&sbuf_), sbuf_(vec) {}
	G3BufferInputStream(char *buf, size_t len) : std::istream(&sbuf_), sbuf_(buf, len) {}
private:
	G3InputStreamBuffer sbuf_;
};

class G3OutputStreamBuffer : public std::basic_streambuf<char>
{
public:
	G3OutputStreamBuffer(std::vector<char> &buffer) : buffer_(buffer) {
		setp(buffer_.data(), buffer_.data() + buffer_.size());
	}
protected:
	std::streamsize xsputn(const char* s, std::streamsize n) {
		buffer_.insert(buffer_.end(), s, s + n);
		pbump(n);
		return n;
	}

	int overflow(int c = std::char_traits<char>::eof()) {
		if (c != std::char_traits<char>::eof()) {
			buffer_.push_back(static_cast<char>(c));
			pbump(1);
		}
		return c;
	}
private:
	std::vector<char>& buffer_;
};

class G3BufferOutputStream : public std::ostream
{
public:
	G3BufferOutputStream(std::vector<char> &vec) : std::ostream(&sbuf_), sbuf_(vec) {}
private:
	G3OutputStreamBuffer sbuf_;
};

// Pickling interface
// For classes with private default constructors, privately add this as a friend class.
class G3Pickler
{
public:
	template <class T> static
	auto dumpstate(const T &obj) {
		std::vector<char> buffer;
		G3BufferOutputStream os(buffer);
		{
			cereal::PortableBinaryOutputArchive ar(os);
			ar << obj;
		}
		os.flush();

		return py::bytes(buffer.data(), buffer.size());
	}

	template <class T> static
	auto getstate(const py::object &self) {
		py::bytes data = dumpstate<T>(self.cast<const T &>());

		py::dict py_state;
		if (py::hasattr(self, "__dict__"))
			py_state = self.attr("__dict__");

		return py::make_tuple(py_state, data);
	}

	template <class T> static
	auto loadstate(T &obj, const std::string_view &buffer) {
		G3BufferInputStream fis((char *)&buffer[0], buffer.size());
		cereal::PortableBinaryInputArchive ar(fis);

		ar >> obj;
	}

	template <class T> static
	auto setstate(const py::tuple &state) {
		auto py_state = state[0].cast<py::dict>();
		auto buffer = state[1].cast<std::string_view>();

		T obj;
		loadstate<T>(obj, buffer);

		return std::make_pair(std::move(obj), py_state);
	}
};

// Call this function in a class .def method to enable pickling
template <class T>
struct g3frameobject_picklesuite
{
	static constexpr bool op_enable_if_hook = true;

	template <typename Class>
	void execute(Class &cl) const {
		cl.def(py::pickle(&G3Pickler::getstate<T>, &G3Pickler::setstate<T>));
		cl.def("_cereal_loads", &G3Pickler::loadstate<T>,
		    "Populate this instance from a serialized data buffer");
		cl.def("_cereal_dumps", &G3Pickler::dumpstate<T>,
		    "Save the state of this instance to a serialized data buffer");
	}
};

// Register a G3FrameObject-derived class.  Includes a copy constructor,
// pickling interface, and string representation.
template <typename T, typename... Bases, typename... Args>
auto
register_frameobject(py::module_ &scope, const std::string &name, Args&&...args)
{
	auto cls = register_class<T, Bases..., G3FrameObject>(scope, name.c_str(),
	    std::forward<Args>(args)...);

	// copy constructor
	cls.def(py::init<const T &>(), "Copy constructor");

	// pickling infrastructure
	cls.def(g3frameobject_picklesuite<T>());

	// string representation
	cls.def("__str__", &T::Summary)
	    .def("Summary", &T::Summary, "Short (one-line) description of the object")
	    .def("Description", &T::Description,
	        "Long-form human-readable description of the object");

	return cls;
}

#endif
