#ifndef _G3_SERIALIZATION_H
#define _G3_SERIALIZATION_H

#include <dataio.h>

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

template <class T>
struct g3frameobject_picklesuite : boost::python::pickle_suite
{
	static boost::python::tuple getstate(boost::python::object obj)
	{
		namespace bp = boost::python;
		std::vector<char> buffer;
		boost::iostreams::filtering_ostream os;
		g3_ostream_to_buffer(os, buffer);
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

		boost::iostreams::filtering_istream fis;
		g3_istream_from_buffer(fis, (char *)view.buf, view.len);
		cereal::PortableBinaryInputArchive ar(fis);

		bp::extract<bp::dict>(obj.attr("__dict__"))().update(state[0]);
		ar >> bp::extract<T &>(obj)();
		PyBuffer_Release(&view);
	}
};

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
