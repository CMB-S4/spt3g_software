#ifndef _G3_SERIALIZATION_H
#define _G3_SERIALIZATION_H

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_stream.hpp>

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
		boost::iostreams::stream<
		    boost::iostreams::back_insert_device<std::vector<char> > >
		    os(buffer);
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

		boost::iostreams::array_source src((char *)view.buf, view.len);
		boost::iostreams::filtering_istream fis(src);
		cereal::PortableBinaryInputArchive ar(fis);

		bp::extract<bp::dict>(obj.attr("__dict__"))().update(state[0]);
		ar >> bp::extract<T &>(obj)();
		PyBuffer_Release(&view);
	}
};

template <class T>
inline int
_g3_class_version(T *)
{
	return cereal::detail::Version<T>::version;
}

#define G3_CHECK_VERSION(v) \
	if (v > _g3_class_version(this)) \
		log_fatal("Trying to read newer class version (%d) than " \
		    "supported (%d). Please upgrade your software.", v, \
		    _g3_class_version(this));
	

#endif
