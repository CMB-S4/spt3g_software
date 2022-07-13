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

G3SkyMapMask::G3SkyMapMask(const G3SkyMapMask &m) : G3FrameObject()
{
	parent_ = m.Parent()->Clone(false);
	data_ = std::vector<bool>(m.data_);
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

G3SkyMapMask &
G3SkyMapMask::operator |=(const G3SkyMapMask &rhs)
{
	g3_assert(IsCompatible(rhs));

	for (size_t i = 0; i < data_.size(); i++)
		(*this)[i] = rhs[i] || (*this)[i];

	return *this;
}

G3SkyMapMask &
G3SkyMapMask::operator &=(const G3SkyMapMask &rhs)
{
	g3_assert(IsCompatible(rhs));

	for (size_t i = 0; i < data_.size(); i++)
		(*this)[i] = rhs[i] && (*this)[i];

	return *this;
}

G3SkyMapMask &
G3SkyMapMask::operator ^=(const G3SkyMapMask &rhs)
{
	g3_assert(IsCompatible(rhs));

	for (size_t i = 0; i < data_.size(); i++)
		(*this)[i] = rhs[i] ^ (*this)[i];

	return *this;
}

G3SkyMapMask &
G3SkyMapMask::invert()
{
	for (size_t i = 0; i < data_.size(); i++)
		(*this)[i] = !(*this)[i];

	return *this;
}

G3SkyMapMask
G3SkyMapMask::operator ~()
{
	G3SkyMapMask mask(*Parent());
	for (size_t i = 0; i < data_.size(); i++)
		mask[i] = !(*this)[i];

	return mask;
}

G3SkyMapMask
G3SkyMapMask::operator |(const G3SkyMapMask &rhs)
{
	g3_assert(IsCompatible(rhs));

	G3SkyMapMask mask(*Parent());
	for (size_t i = 0; i < data_.size(); i++)
		mask[i] = (*this)[i] || rhs[i];

	return mask;
}

G3SkyMapMask
G3SkyMapMask::operator &(const G3SkyMapMask &rhs)
{
	g3_assert(IsCompatible(rhs));

	G3SkyMapMask mask(*Parent());
	for (size_t i = 0; i < data_.size(); i++)
		mask[i] = (*this)[i] && rhs[i];

	return mask;
}

G3SkyMapMask
G3SkyMapMask::operator ^(const G3SkyMapMask &rhs)
{
	g3_assert(IsCompatible(rhs));

	G3SkyMapMask mask(*Parent());
	for (size_t i = 0; i < data_.size(); i++)
		mask[i] = (*this)[i] ^ rhs[i];

	return mask;
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

static G3SkyMapMaskPtr
skymapmask_pyinvert(G3SkyMapMaskPtr m)
{
	// Reference-counting problems returning references to Python,
	// so use shared pointers.
	m->invert();
	return m;
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
	  .def("invert", &skymapmask_pyinvert, "Invert all elements in mask")
	  .def(bp::self |= bp::self)
	  .def(bp::self &= bp::self)
	  .def(bp::self ^= bp::self)
	  .def(~bp::self)
	  .def(bp::self | bp::self)
	  .def(bp::self & bp::self)
	  .def(bp::self ^ bp::self)
	;
}

