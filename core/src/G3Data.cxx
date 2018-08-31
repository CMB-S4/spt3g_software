#include <pybindings.h>
#include <serialization.h>
#include <G3Data.h>
#include <sstream>

std::string G3Bool::Description() const
{
	if (value)
		return "True";
	else
		return "False";
}

std::string G3Int::Description() const
{
	std::ostringstream s;
	s << value;
	return s.str();
}

std::string G3Double::Description() const
{
	std::ostringstream s;
	s << value;
	return s.str();
}

std::string G3String::Description() const
{
	std::ostringstream s;
	s << "\"" << value << "\"";
	return s.str();
}

template <class A> void G3Bool::serialize(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("value", value);
}

template <class A> void G3Int::serialize(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("value", value);
}

template <class A> void G3Double::serialize(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("value", value);
}

template <class A> void G3String::serialize(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("value", value);
}

G3_SERIALIZABLE_CODE(G3Bool);
G3_SERIALIZABLE_CODE(G3Int);
G3_SERIALIZABLE_CODE(G3Double);
G3_SERIALIZABLE_CODE(G3String);

PYBINDINGS("core") {
	EXPORT_FRAMEOBJECT(G3Bool, init<bool>(), "Serializable boolean type")
	    .def_readwrite("value", &G3Bool::value)
	    .def("__nonzero__", &G3Bool::truth)
	    .def("__bool__", &G3Bool::truth)
	;

	EXPORT_FRAMEOBJECT(G3Int, init<int64_t>(), "Serializable integer type")
	    .def_readwrite("value", &G3Int::value)
	;

	EXPORT_FRAMEOBJECT(G3Double, init<double>(), "Serializable double")
	    .def_readwrite("value", &G3Double::value)
	;

	EXPORT_FRAMEOBJECT(G3String, init<std::string>(), "Serializable string")
	    .def_readwrite("value", &G3String::value)
	;
}

