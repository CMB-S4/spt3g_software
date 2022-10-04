#include <pybindings.h>
#include <serialization.h>
#include <typeinfo>
#include <sys/types.h>

#include <maps/G3SkyMapMask.h>
#include <maps/G3SkyMap.h>
#include <maps/FlatSkyMap.h>

G3SkyMapMask::G3SkyMapMask(const G3SkyMap &parent, bool use_data,
  bool zero_nans, bool zero_infs) : G3FrameObject()
{
	// masks have no units
	G3SkyMapPtr tmp = parent.Clone(false);
	tmp->units = G3Timestream::None;
	tmp->pol_type = G3SkyMap::None;
	tmp->SetPolConv(G3SkyMap::ConvNone);
	tmp->weighted = false;
	parent_ = tmp;
	data_ = std::vector<bool>(parent.size());

	if (use_data)
		FillFromMap(parent, zero_nans, zero_infs);
}

G3SkyMapMask::G3SkyMapMask(const G3SkyMap &parent, boost::python::object v,
  bool zero_nans, bool zero_infs) : G3FrameObject()
{
	// masks have no units
	G3SkyMapPtr tmp = parent.Clone(false);
	tmp->units = G3Timestream::None;
	tmp->pol_type = G3SkyMap::None;
	tmp->SetPolConv(G3SkyMap::ConvNone);
	tmp->weighted = false;
	parent_ = tmp;
	data_ = std::vector<bool>(parent.size());

	// fall back to simple constructor
	if (boost::python::extract<bool>(v).check()) {
		if (boost::python::extract<bool>(v)())
			FillFromMap(parent, zero_nans, zero_infs);
		return;
	}

	FillFromArray(v, zero_nans, zero_infs);
}

G3SkyMapMask::G3SkyMapMask(const G3SkyMapMask &m) : G3FrameObject()
{
	parent_ = m.parent_;
	data_ = std::vector<bool>(m.data_);
}

G3SkyMapMaskPtr
G3SkyMapMask::Clone(bool copy_data) const
{
	if (copy_data)
		return boost::make_shared<G3SkyMapMask>(*this);
	else
		return boost::make_shared<G3SkyMapMask>(*Parent());
}

void
G3SkyMapMask::FillFromMap(const G3SkyMap &m, bool zero_nans, bool zero_infs)
{
	g3_assert(IsCompatible(m));

	for (size_t i = 0; i < m.size(); i++) {
		double v = m.at(i);
		if (v == 0)
			continue;
		if (zero_nans && std::isnan(v))
			continue;
		if (zero_infs && std::isinf(v))
			continue;
		data_[i] = true;
	}
}

#define FILL_BUFFER(type) \
	for (size_t i = 0; i < view.len / sizeof(type); i++) { \
		double v = ((type *)view.buf)[i]; \
		if (v == 0) \
			continue; \
		if (zero_nans && std::isnan(v)) \
			continue; \
		if (zero_infs && std::isinf(v)) \
			continue; \
		data_[i] = true; \
	}

void
G3SkyMapMask::FillFromArray(boost::python::object v, bool zero_nans, bool zero_infs)
{

	Py_buffer view;

	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_C_CONTIGUOUS) != -1) {
		if (view.ndim != 1) {
			PyBuffer_Release(&view);
			log_fatal("Only 1-D masks supported");
		}

		size_t npix = view.shape[0];
		if (npix != size()) {
			PyBuffer_Release(&view);
			log_fatal("Got array of shape (%zu,), expected (%zu,)", npix, size());
		}

		const char *format = view.format;

		if (format[0] == '@' || format[0] == '=')
			format++;
#if BYTE_ORDER == LITTLE_ENDIAN
		else if (format[0] == '<')
			format++;
		else if (format[0] == '>' || format[0] == '!') {
			PyBuffer_Release(&view);
			log_fatal("Does not support big-endian numpy arrays");
		}
#else
		else if (format[0] == '<') {
			PyBuffer_Release(&view);
			log_fatal("Does not support little-endian numpy arrays");
		} else if (format[0] == '>' || format[0] == '!')
			format++;
#endif

		if (strcmp(format, "d") == 0) {
			FILL_BUFFER(double);
		} else if (strcmp(format, "f") == 0) {
			FILL_BUFFER(float);
		} else if (strcmp(format, "i") == 0) {
			FILL_BUFFER(int);
		} else if (strcmp(format, "I") == 0) {
			FILL_BUFFER(unsigned int);
		} else if (strcmp(format, "l") == 0) {
			FILL_BUFFER(long);
		} else if (strcmp(format, "L") == 0) {
			FILL_BUFFER(unsigned long);
		} else if (strcmp(format, "b") == 0) {
			FILL_BUFFER(char);
		} else if (strcmp(format, "B") == 0) {
			FILL_BUFFER(unsigned char);
		} else if (strcmp(format, "?") == 0) {
			FILL_BUFFER(bool);
		} else {
			PyBuffer_Release(&view);
			log_fatal("Unknown type code %s", view.format);
		}
		PyBuffer_Release(&view);

		return;
	}

	throw boost::python::error_already_set();
}

G3SkyMapMaskPtr
G3SkyMapMask::ArrayClone(boost::python::object v, bool zero_nans, bool zero_infs) const
{
	G3SkyMapMaskPtr m = Clone(false);
	m->FillFromArray(v, zero_nans, zero_infs);
	return m;
}

G3SkyMapMask::const_iterator::const_iterator(const G3SkyMapMask &mask, bool begin) :
    mask_(mask)
{
	index_ = begin ? 0 : mask_.size();
	set_value();
}

G3SkyMapMask::const_iterator
G3SkyMapMask::const_iterator::operator++()
{
	++index_;
	set_value();

	return *this;
}

size_t
G3SkyMapMask::size() const
{
	return data_.size();
}

std::vector<bool>::reference
G3SkyMapMask::operator [] (size_t i)
{
	return data_[i];
}

bool
G3SkyMapMask::at (size_t i) const
{
	if (i >= data_.size())
		return false;
	return data_[i];
}

G3SkyMapMask &
G3SkyMapMask::operator |=(const G3SkyMapMask &rhs)
{
	g3_assert(IsCompatible(rhs));

	for (size_t i = 0; i < size(); i++)
		(*this)[i] = rhs.at(i) || at(i);

	return *this;
}

G3SkyMapMask &
G3SkyMapMask::operator &=(const G3SkyMapMask &rhs)
{
	g3_assert(IsCompatible(rhs));

	for (auto i: *this)
		(*this)[i.first] = rhs.at(i.first) && i.second;

	return *this;
}

G3SkyMapMask &
G3SkyMapMask::operator ^=(const G3SkyMapMask &rhs)
{
	g3_assert(IsCompatible(rhs));

	for (size_t i = 0; i < size(); i++)
		(*this)[i] = rhs.at(i) ^ at(i);

	return *this;
}

G3SkyMapMask &
G3SkyMapMask::invert()
{
	for (size_t i = 0; i < size(); i++)
		(*this)[i] = !at(i);

	return *this;
}

bool
G3SkyMapMask::all() const
{
	for (size_t i = 0; i < size(); i++)
		if (at(i) == 0)
			return false;
	return true;
}

bool
G3SkyMapMask::any() const
{
	for (auto i: *this)
		if (i.second != 0)
			return true;
	return false;
}

size_t
G3SkyMapMask::sum() const
{
	size_t s = 0;
	for (auto i: *this)
		s += i.second;
	return s;
}

std::vector<uint64_t>
G3SkyMapMask::NonZeroPixels() const
{
	std::vector<uint64_t> indices;

	for (auto i: *this) {
		if (i.second)
			indices.push_back(i.first);
	}

	return indices;
}

G3SkyMapMask
G3SkyMapMask::operator ~() const
{
	G3SkyMapMask mask(*Parent());
	for (size_t i = 0; i < size(); i++)
		if (!at(i))
			mask[i] = true;

	return mask;
}

G3SkyMapMask
G3SkyMapMask::operator |(const G3SkyMapMask &rhs) const
{
	g3_assert(IsCompatible(rhs));

	G3SkyMapMask mask(*Parent());
	for (size_t i = 0; i < size(); i++)
		if (at(i) || rhs.at(i))
			mask[i] = true;

	return mask;
}

G3SkyMapMask
G3SkyMapMask::operator &(const G3SkyMapMask &rhs) const
{
	g3_assert(IsCompatible(rhs));

	G3SkyMapMask mask(*Parent());
	for (auto i: *this)
		if (i.second && rhs.at(i.first))
			mask[i.first] = true;

	return mask;
}

G3SkyMapMask
G3SkyMapMask::operator ^(const G3SkyMapMask &rhs) const
{
	g3_assert(IsCompatible(rhs));

	G3SkyMapMask mask(*Parent());
	for (size_t i = 0; i < size(); i++)
		if (at(i) ^ rhs.at(i))
			mask[i] = true;

	return mask;
}

G3SkyMapMask
G3SkyMapMask::operator ==(const G3SkyMapMask &rhs) const
{
	g3_assert(IsCompatible(rhs));

	G3SkyMapMask mask(*Parent());
	for (size_t i = 0; i < size(); i++)
		if (at(i) == rhs.at(i))
			mask[i] = true;

	return mask;
}

G3SkyMapMask
G3SkyMapMask::operator !=(const G3SkyMapMask &rhs) const
{
	g3_assert(IsCompatible(rhs));

	G3SkyMapMask mask(*Parent());
	for (size_t i = 0; i < size(); i++)
		if (at(i) != rhs.at(i))
			mask[i] = true;

	return mask;
}

void
G3SkyMapMask::ApplyMask(const G3SkyMapMask &rhs, bool inverse)
{
	g3_assert(IsCompatible(rhs));

	for (auto i: *this)
		if (i.second && rhs.at(i.first) == inverse)
			(*this)[i.first] = false;
}

bool
G3SkyMapMask::IsCompatible(const G3SkyMap &map) const
{
	return Parent()->IsCompatible(map);
}

bool
G3SkyMapMask::IsCompatible(const G3SkyMapMask &mask) const
{
	return Parent()->IsCompatible(*mask.Parent());
}

G3SkyMapPtr
G3SkyMapMask::MakeBinaryMap() const
{
	G3SkyMapPtr out = Parent()->Clone();

	for (auto i: *this) {
		if (i.second)
			(*out)[i.first] = 1.0;
	}

	return out;
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
skymapmask_getitem(const G3SkyMapMask &m, boost::python::object index)
{
	using namespace boost::python;

	int i = 0;
	if (extract<int>(index).check()) {
		i = extract<int>(index)();

		if (i < 0)
			i = m.size() + i;
	} else if (extract<tuple>(index).check()) {
		int x, y;

		tuple t = extract<tuple>(index)();
		FlatSkyMapConstPtr fsm = boost::dynamic_pointer_cast<const FlatSkyMap>(m.Parent());
		if (!fsm) {
			PyErr_SetString(PyExc_TypeError,
			    "N-D pixels, but underlying map is not a flat sky map");
			boost::python::throw_error_already_set();
		}

		x = extract<int>(t[1])();
		y = extract<int>(t[0])();
		if (x < 0)
			x += fsm->shape()[0];
		if (y < 0)
			y += fsm->shape()[0];
		if (size_t(x) >= fsm->shape()[0] ||
		    size_t(y) >= fsm->shape()[1]) { 
			PyErr_SetString(PyExc_IndexError, "Index out of range");
			boost::python::throw_error_already_set();
		}

		i = y * fsm->shape()[0] + x;
	} else {
		PyErr_SetString(PyExc_TypeError,
		    "Need to pass an integer pixel ID or (optionally) for 2D maps a tuple of coordinates");
		boost::python::throw_error_already_set();
	}
	
	if (i < 0 || size_t(i) >= m.size()) {
		PyErr_SetString(PyExc_IndexError, "Index out of range");
		boost::python::throw_error_already_set();
	}

	return m.at(i);
}

static void
skymapmask_setitem(G3SkyMapMask &m, boost::python::object index, bool val)
{
	using namespace boost::python;

	int i = 0;
	if (extract<int>(index).check()) {
		i = extract<int>(index)();

		if (i < 0)
			i = m.size() + i;
	} else if (extract<tuple>(index).check()) {
		int x, y;

		tuple t = extract<tuple>(index)();
		FlatSkyMapConstPtr fsm = boost::dynamic_pointer_cast<const FlatSkyMap>(m.Parent());
		if (!fsm) {
			PyErr_SetString(PyExc_TypeError,
			    "N-D pixels, but underlying map is not a flat sky map");
			boost::python::throw_error_already_set();
		}

		x = extract<int>(t[1])();
		y = extract<int>(t[0])();
		if (x < 0)
			x += fsm->shape()[0];
		if (y < 0)
			y += fsm->shape()[0];
		if (size_t(x) >= fsm->shape()[0] ||
		    size_t(y) >= fsm->shape()[1]) { 
			PyErr_SetString(PyExc_IndexError, "Index out of range");
			boost::python::throw_error_already_set();
		}

		i = y * fsm->shape()[0] + x;
	} else {
		PyErr_SetString(PyExc_TypeError,
		    "Need to pass an integer pixel ID or (optionally) for 2D maps a tuple of coordinates");
		boost::python::throw_error_already_set();
	}

	if (i < 0 || size_t(i) >= m.size()) {
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

// This function handles some implicit pointer conversions that boost won't
// do because G3SkyMap is an abstract type. Just do it by hand.
static G3SkyMapPtr
skymapmask_pyparent(G3SkyMapMaskPtr m)
{
	return boost::const_pointer_cast<G3SkyMap>(m->Parent());
}

static bool
skymapmask_pybool(G3SkyMapMaskPtr m)
{
	PyErr_SetString(PyExc_ValueError,
	    "ValueError: The truth value of a G3SkyMapMask is ambiguous. Use m.any() or m.all()");
	bp::throw_error_already_set();

	return false;
}

PYBINDINGS("maps")
{
	using namespace boost::python;

	EXPORT_FRAMEOBJECT(G3SkyMapMask, no_init,
	    "Boolean mask of a sky map. Set pixels to use to true, pixels to "
	    "ignore to false. If use_data set in contrast, mask initialized to "
	    "true where input map is non-zero; otherwise, all elements are "
	    "initialized to zero.  Use zero_nans or zero_infs to exclude nan "
	    "or inf elements from the mask.")
	  .def(bp::init<const G3SkyMap &, bool, bool, bool>(
	       (bp::arg("parent"),
		bp::arg("use_data")=false,
		bp::arg("zero_nans")=false,
		bp::arg("zero_infs")=false),
	    "Instantiate a G3SkyMapMask from a parent G3SkyMap"))
	  .def(bp::init<const G3SkyMap &, bp::object, bool, bool>(
	       (bp::arg("parent"),
		bp::arg("data"),
		bp::arg("zero_nans")=false,
		bp::arg("zero_infs")=false),
	    "Instantiate a G3SkyMapMask from a 1D numpy array"))
	  .def("clone", &G3SkyMapMask::Clone, (bp::arg("copy_data")=true),
	    "Return a mask of the same type, populated with a copy of the data "
	    "if the argument is true (default), empty otherwise.")
	  .def("array_clone", &G3SkyMapMask::ArrayClone,
	    (bp::arg("data"), bp::arg("zero_nans")=false, bp::arg("zero_infs")=false),
	    "Return a mask of the same type, populated from the input numpy array")
	  .add_property("parent", &skymapmask_pyparent, "\"Parent\" map which "
	    "contains no data, but can be used to retrieve the parameters of "
	    "the map to which this mask corresponds.")
	  .add_property("size", &G3SkyMapMask::size, "Number of pixels in mask")
	  .def<bool (G3SkyMapMask::*)(const G3SkyMapMask &) const>("compatible", &G3SkyMapMask::IsCompatible, "Returns true if the two masks can be applied to the same map.")
	  .def<bool (G3SkyMapMask::*)(const G3SkyMap &) const>("compatible", &G3SkyMapMask::IsCompatible, "Returns true if this mask can be applied to the given map.")
	  .def("__getitem__", &skymapmask_getitem)
	  .def("__setitem__", &skymapmask_setitem)
	  .def("__bool__", &skymapmask_pybool)
	  .def("invert", &skymapmask_pyinvert, "Invert all elements in mask")
	  .def("_call", &G3SkyMapMask::all, "Test whether all elements are non-zero")
	  .def("_cany", &G3SkyMapMask::any, "Test whether any elements are non-zero")
	  .def("_csum", &G3SkyMapMask::sum, "Sum of all elements in mask")
	  .def("nonzero", &G3SkyMapMask::NonZeroPixels,
	       "Return a list of indices of non-zero pixels in the mask")
	  .def("apply_mask", &G3SkyMapMask::ApplyMask,
	    (bp::arg("mask"), bp::arg("inverse")=false),
	    "Apply a mask in-place to the mask, optionally inverting which "
	    "pixels are zeroed.  If inverse = False, this is equivalent to "
	    "in-place element-wise logical-and with the mask.")
	  .def(bp::self |= bp::self)
	  .def(bp::self &= bp::self)
	  .def(bp::self ^= bp::self)
	  .def(~bp::self)
	  .def(bp::self | bp::self)
	  .def(bp::self & bp::self)
	  .def(bp::self ^ bp::self)
	  .def(bp::self == bp::self)
	  .def(bp::self != bp::self)
	  .def("to_map", &G3SkyMapMask::MakeBinaryMap, "Create a skymap with data set to the contents of this mask (1.0 where True, 0.0 where False), which can be useful for plotting.")
	;
	register_pointer_conversions<G3SkyMapMask>();
}

