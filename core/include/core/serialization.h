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
		setg(&vec[0], &vec[0], &vec[0] + vec.size());
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
		setp(&buffer_[0], &buffer_[0] + buffer_.size());
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

template <class T>
struct g3frameobject_picklesuite : py::pickle_suite
{
	static py::tuple getstate(const py::object &obj)
	{
		std::vector<char> buffer;
		G3BufferOutputStream os(buffer);
		{
			cereal::PortableBinaryOutputArchive ar(os);
			ar << py::extract<const T &>(obj)();
		}
		os.flush();

		return py::make_tuple(obj.attr("__dict__"),
		    py::object(py::handle<>(
		    PyBytes_FromStringAndSize(&buffer[0], buffer.size()))));
	}

	static void setstate(py::object &obj, py::tuple state)
	{
		Py_buffer view;
		PyObject_GetBuffer(py::object(state[1]).ptr(), &view,
		    PyBUF_SIMPLE);

		G3BufferInputStream fis((char *)view.buf, view.len);
		cereal::PortableBinaryInputArchive ar(fis);

		py::extract<py::dict>(obj.attr("__dict__"))().update(state[0]);
		ar >> py::extract<T &>(obj)();
		PyBuffer_Release(&view);
	}
};

template <typename T>
void
register_pointer_conversions()
{
	py::implicitly_convertible<std::shared_ptr<T>, G3FrameObjectPtr>();
	py::implicitly_convertible<std::shared_ptr<T>, std::shared_ptr<const T> >();
	py::implicitly_convertible<std::shared_ptr<T>, G3FrameObjectConstPtr>();
}

// Register a G3FrameObject-derived class.  Includes a copy constructor,
// pickling interface, and string representation.
template <typename T, typename... Bases, typename... Args>
auto
register_frameobject(py::module_ &scope, const std::string &name, Args&&...args)
{
	auto cls = register_class<T, Bases..., G3FrameObject>(scope, name.c_str(),
	    std::forward<Args>(args)...);

	// copy constructor
	cls.def(py::init<const T &>("Copy constructor"));

	// pickling infrastructure
	cls.def_pickle(g3frameobject_picklesuite<T>());

	// string representation
	cls.def("__str__", &T::Summary)
	    .def("Summary", &T::Summary, "Short (one-line) description of the object")
	    .def("Description", &T::Description,
	        "Long-form human-readable description of the object");

	register_pointer_conversions<T>();

	return cls;
}

#endif
