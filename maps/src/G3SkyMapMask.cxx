#include <pybindings.h>
#include <serialization.h>
#include <container_pybindings.h>
#include <typeinfo>
#include <algorithm>
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
	tmp->pol_conv = G3SkyMap::ConvNone;
	tmp->weighted = false;
	parent_ = tmp;
	data_ = std::vector<bool>(parent.size());

	if (use_data) {
		g3_assert(IsCompatible(parent));

		for (size_t i = 0; i < parent.size(); i++) {
			double v = parent.at(i);
			if (v == 0)
				continue;
			if (zero_nans && std::isnan(v))
				continue;
			if (zero_infs && std::isinf(v))
				continue;
			data_[i] = true;
		}
	}
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
		return std::make_shared<G3SkyMapMask>(*this);
	else
		return std::make_shared<G3SkyMapMask>(*Parent());
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
void G3SkyMapMask::load(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("parent", parent_);
	if (v < 2) {
		ar & make_nvp("data", data_);
	} else {
		std::vector<uint8_t> packed;
		size_t nbits;
		ar & make_nvp("data", packed);
		data_.resize(packed.size()*8);
		for (size_t i = 0; i < packed.size(); i++)
			for (int j = 0; j < 8; j++)
				data_[i*8 + j] = (packed[i] >> j) & 1;
		ar & make_nvp("nbits", nbits); // In case not a multiple of 8
		data_.resize(nbits);
	}
}

template <class A>
void G3SkyMapMask::save(A &ar, unsigned v) const
{
	using namespace cereal;
	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("parent", parent_);

	std::vector<uint8_t> packed((data_.size() / 8) +
	    ((data_.size() % 8) ? 1 : 0));
	// Pack data in all bits up to a multiple of 8, then the remainder
	// Two pieces so the compiler can unroll the inner loop in the first
	// part
	for (size_t i = 0; i < data_.size() / 8; i++) {
		packed[i] = 0;
		for (int j = 0; j < 8; j++)
			packed[i] |= !!data_[i*8 + j] << j;
	}
	if (data_.size() % 8) {
		const ssize_t i = packed.size() - 1;
		packed[i] = 0;
		for (int j = 0; j < (ssize_t) data_.size() - i*8; j++)
			packed[i] |= !!data_[i*8 + j] << j;
	}

	ar & make_nvp("data", packed);
	ar & make_nvp("nbits", data_.size());
}
		

G3_SPLIT_SERIALIZABLE_CODE(G3SkyMapMask);

#define FILL_BUFFER(type) \
	for (size_t i = 0; i < view.len / sizeof(type); i++) { \
		double v = ((type *)view.buf)[i]; \
		if (v == 0) \
			continue; \
		if (zero_nans && std::isnan(v)) \
			continue; \
		if (zero_infs && std::isinf(v)) \
			continue; \
		(*mask)[i] = true; \
	}

static G3SkyMapMaskPtr
skymapmask_from_numpy(const G3SkyMap &parent, py::object v,
  bool zero_nans, bool zero_infs)
{
	// fall back to simple constructor
	auto ext  = py::extract<bool>(v);
	if (ext.check())
		return std::make_shared<G3SkyMapMask>(parent, ext(), zero_nans, zero_infs);

	G3SkyMapMaskPtr mask(new G3SkyMapMask(parent));

	Py_buffer view;

	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_C_CONTIGUOUS) == -1)
		throw py::error_already_set();

	if (view.ndim != 1) {
		PyBuffer_Release(&view);
		log_fatal("Only 1-D masks supported");
	}

	size_t npix = view.shape[0];
	if (npix != mask->size()) {
		PyBuffer_Release(&view);
		log_fatal("Got array of shape (%zu,), expected (%zu,)",
		    npix, mask->size());
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

	return mask;
}

static G3SkyMapMaskPtr
skymapmask_array_clone(const G3SkyMapMask &m, py::object v,
    bool zero_nans, bool zero_infs)
{
	return skymapmask_from_numpy(*m.Parent(), v, zero_nans, zero_infs);
}

static int
skymapmask_index(const G3SkyMapMask &m, const py::object &index)
{
	int i = 0;
	if (py::extract<int>(index).check()) {
		i = py::extract<int>(index)();

		if (i < 0)
			i = m.size() + i;
	} else if (py::extract<py::tuple>(index).check()) {
		int x, y;

		py::tuple t = py::extract<py::tuple>(index)();
		FlatSkyMapConstPtr fsm = std::dynamic_pointer_cast<const FlatSkyMap>(m.Parent());
		if (!fsm)
			throw py::type_error("N-D pixels, but underlying map is not a flat sky map");

		x = py::extract<int>(t[1])();
		y = py::extract<int>(t[0])();
		if (x < 0)
			x += fsm->shape()[0];
		if (y < 0)
			y += fsm->shape()[0];
		if (size_t(x) >= fsm->shape()[0] ||
		    size_t(y) >= fsm->shape()[1])
			throw py::index_error("Index out of range");

		i = y * fsm->shape()[0] + x;
	} else {
		throw py::type_error("Need to pass an integer pixel ID or "
		    "(optionally) for 2D maps a tuple of coordinates");
	}
	
	if (i < 0 || size_t(i) >= m.size())
		throw py::index_error("Index out of range");

	return i;
}

static bool
skymapmask_getitem(const G3SkyMapMask &m, py::object index)
{
	return m.at(skymapmask_index(m, index));
}

static void
skymapmask_setitem(G3SkyMapMask &m, py::object index, bool val)
{
	m[skymapmask_index(m, index)] = val;
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
	return std::const_pointer_cast<G3SkyMap>(m->Parent());
}

static bool
skymapmask_pybool(G3SkyMapMaskPtr m)
{
	throw py::value_error("The truth value of a G3SkyMapMask "
	    "is ambiguous. Use m.any() or m.all()");

	return false;
}

static py::dict
G3SkyMapMask_array_interface(const G3SkyMapMask &self)
{
	py::dict out;
	out["typestr"] = "f8";
	out["data"] = py::object(self.MakeBinaryMap());
	auto shape = self.Parent()->shape();
	std::reverse(shape.begin(), shape.end());
	std::vector<uint64_t> ushape(shape.begin(), shape.end());
	out["shape"] = py::tuple(ushape);

	return out;
}

PYBINDINGS("maps", scope)
{
	EXPORT_FRAMEOBJECT(G3SkyMapMask, py::no_init,
	    "Boolean mask of a sky map. Set pixels to use to true, pixels to "
	    "ignore to false. If use_data set in contrast, mask initialized to "
	    "true where input map is non-zero; otherwise, all elements are "
	    "initialized to zero.  Use zero_nans or zero_infs to exclude nan "
	    "or inf elements from the mask.")
	  .def(py::init<const G3SkyMap &, bool, bool, bool>(
	       (py::arg("parent"),
	       py::arg("use_data")=false,
	       py::arg("zero_nans")=false,
	       py::arg("zero_infs")=false),
	    "Instantiate a G3SkyMapMask from a parent G3SkyMap"))
	  .def("__init__", py::make_constructor(skymapmask_from_numpy, py::default_call_policies(),
	       (py::arg("parent"),
	       py::arg("data"),
	       py::arg("zero_nans")=false,
	       py::arg("zero_infs")=false)),
	    "Instantiate a G3SkyMapMask from a 1D numpy array")
	  .def("clone", &G3SkyMapMask::Clone, (py::arg("copy_data")=true),
	    "Return a mask of the same type, populated with a copy of the data "
	    "if the argument is true (default), empty otherwise.")
	  .def("array_clone", &skymapmask_array_clone,
	    (py::arg("data"), py::arg("zero_nans")=false, py::arg("zero_infs")=false),
	    "Return a mask of the same type, populated from the input numpy array")
	  .add_property("parent", &skymapmask_pyparent, "\"Parent\" map which "
	    "contains no data, but can be used to retrieve the parameters of "
	    "the map to which this mask corresponds.")
	  .add_property("size", &G3SkyMapMask::size, "Number of pixels in mask")
	  .def("compatible", (bool (G3SkyMapMask::*)(const G3SkyMapMask &) const)
	    &G3SkyMapMask::IsCompatible,
	    "Returns true if the two masks can be applied to the same map.")
	  .def("compatible", (bool (G3SkyMapMask::*)(const G3SkyMap &) const)
	    &G3SkyMapMask::IsCompatible,
	    "Returns true if this mask can be applied to the given map.")
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
	    (py::arg("mask"), py::arg("inverse")=false),
	    "Apply a mask in-place to the mask, optionally inverting which "
	    "pixels are zeroed.  If inverse = False, this is equivalent to "
	    "in-place element-wise logical-and with the mask.")
	  .def(py::self |= py::self)
	  .def(py::self &= py::self)
	  .def(py::self ^= py::self)
	  .def(~py::self)
	  .def(py::self | py::self)
	  .def(py::self & py::self)
	  .def(py::self ^ py::self)
	  .def(py::self == py::self)
	  .def(py::self != py::self)
	  .add_property("__array_interface__", G3SkyMapMask_array_interface)
	  .def("to_map", &G3SkyMapMask::MakeBinaryMap,
	    "Create a skymap with data set to the contents of this mask "
	    "(1.0 where True, 0.0 where False), which can be useful for plotting.")
	;
	register_pointer_conversions<G3SkyMapMask>();
}

