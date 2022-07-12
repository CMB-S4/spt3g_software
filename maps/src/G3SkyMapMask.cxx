#include <pybindings.h>
#include <serialization.h>
#include <typeinfo>
#include <sys/types.h>

#include <maps/G3SkyMapMask.h>

G3SkyMapMask::G3SkyMapMask(const G3SkyMap &parent) : G3FrameObject()
{
	parent_ = parent.Clone(false);
	data_ = std::vector<bool>(parent.size());
}

std::vector<bool>::reference
G3SkyMapMask::operator [] (size_t i)
{
	return data_[i];
}

bool
G3SkyMapMask::at (size_t i) const
{
	return data_[i];
}

template <class A>
void G3SkyMapMask::serialize(A &ar, unsigned v)
{
	using namespace cereal;
	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("parent", parent_);
	ar & make_nvp("data", data_);
}

G3_SERIALIZABLE_CODE(G3SkyMapMask);

static bool
skymapmask_getitem(const G3SkyMapMask &m, int i)
{
	if (i < 0)
		i = m.Parent()->size() + i;
	if (size_t(i) >= m.Parent()->size()) {
		PyErr_SetString(PyExc_IndexError, "Index out of range");
		boost::python::throw_error_already_set();
	}

	return m.at(i);
}

static void
skymapmask_setitem(G3SkyMapMask &m, int i, bool val)
{
	if (i < 0)
		i = m.Parent()->size() + i;
	if (size_t(i) >= m.Parent()->size()) {
		PyErr_SetString(PyExc_IndexError, "Index out of range");
		boost::python::throw_error_already_set();
	}

	m[i] = val;
}

PYBINDINGS("maps")
{
	using namespace boost::python;

	EXPORT_FRAMEOBJECT(G3SkyMapMask, init<const G3SkyMap &>(),
	    "Boolean mask of a sky map. Set pixels to use to true, pixels to "
	    "ignore to false.")
	  .def_readonly("parent", &G3SkyMapMask::Parent, "\"Parent\" map which "
	    "contains no data, but can be used to retrieve the parameters of "
	    "the map to which this mask corresponds.")
	  .def("__getitem__", &skymapmask_getitem)
	  .def("__setitem__", &skymapmask_setitem)
	;
}

