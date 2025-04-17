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
struct g3frameobject_picklesuite : boost::python::pickle_suite
{
	static boost::python::tuple getstate(boost::python::object obj)
	{
		namespace bp = boost::python;
		std::vector<char> buffer;
		G3BufferOutputStream os(buffer);
		{
			cereal::PortableBinaryOutputArchive ar(os);
			ar << bp::extract<const T &>(obj)();
		}
		os.flush();

		return boost::python::make_tuple(obj.attr("__dict__"),
		    bp::object(bp::handle<>(
		    PyBytes_FromStringAndSize(&buffer[0], buffer.size()))));
	}

	static void setstate(boost::python::object obj,	
	    boost::python::tuple state)
	{
		namespace bp = boost::python;
		Py_buffer view;
		PyObject_GetBuffer(bp::object(state[1]).ptr(), &view,
		    PyBUF_SIMPLE);

		G3BufferInputStream fis((char *)view.buf, view.len);
		cereal::PortableBinaryInputArchive ar(fis);

		bp::extract<bp::dict>(obj.attr("__dict__"))().update(state[0]);
		ar >> bp::extract<T &>(obj)();
		PyBuffer_Release(&view);
	}
};

#define EXPORT_FRAMEOBJECT(T, initf, docstring) \
	boost::python::class_<T, boost::python::bases<G3FrameObject>, std::shared_ptr<T> >(#T, docstring, boost::python::initf) \
	    .def(boost::python::init<const T &>()) \
	    .def_pickle(g3frameobject_picklesuite<T>())

#define EXPORT_FRAMEOBJECT_NOINITNAMESPACE(T, initf, docstring) \
	boost::python::class_<T, boost::python::bases<G3FrameObject>, std::shared_ptr<T> >(#T, docstring, initf) \
	    .def(boost::python::init<const T &>()) \
	    .def_pickle(g3frameobject_picklesuite<T>())

template <class T>
inline uint32_t
_g3_class_version(T *)
{
	return cereal::detail::Version<T>::version;
}

#define G3_CHECK_VERSION(v) \
	if ((uint32_t)v > _g3_class_version(this)) \
		log_fatal("Trying to read newer class version (%d) than " \
		    "supported (%d). Please upgrade your software.", v, \
		    _g3_class_version(this));
	

#endif
