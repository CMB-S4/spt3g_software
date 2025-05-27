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

PYBINDINGS("core", scope) {
	register_frameobject<G3Bool>(scope, "G3Bool", "Serializable boolean type")
	    .def(py::init<bool>())
	    .def_readwrite("value", &G3Bool::value)
	    .def("__nonzero__", &G3Bool::truth)
	    .def("__bool__", &G3Bool::truth)
	    .def(py::self == py::self)
	    .def(py::self != py::self)
	;

	register_frameobject<G3Int>(scope, "G3Int", "Serializable integer type")
	    .def(py::init<int64_t>())
	    .def_readwrite("value", &G3Int::value)
	    .def(py::self == py::self)
	    .def(py::self != py::self)
	    .def(py::self >= py::self)
	    .def(py::self > py::self)
	    .def(py::self <= py::self)
	    .def(py::self < py::self)
	;

	register_frameobject<G3Double>(scope, "G3Double", "Serializable double")
	    .def(py::init<double>())
	    .def_readwrite("value", &G3Double::value)
	    .def(py::self == py::self)
	    .def(py::self != py::self)
	    .def(py::self >= py::self)
	    .def(py::self > py::self)
	    .def(py::self <= py::self)
	    .def(py::self < py::self)
	;

	register_frameobject<G3String>(scope, "G3String", "Serializable string")
	    .def(py::init<std::string>())
	    .def_readwrite("value", &G3String::value)
	    .def(py::self == py::self)
	    .def(py::self != py::self)
	    .def(py::self >= py::self)
	    .def(py::self > py::self)
	    .def(py::self <= py::self)
	    .def(py::self < py::self)
	;
}

