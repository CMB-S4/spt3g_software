#include <pybindings.h>
#include <serialization.h>
#include <container_pybindings.h>
#include <G3Vector.h>
#include <complex>

G3_SPLIT_SERIALIZABLE_CODE(G3VectorInt);
G3_SERIALIZABLE_CODE(G3VectorBool);
G3_SERIALIZABLE_CODE(G3VectorDouble);
G3_SERIALIZABLE_CODE(G3VectorComplexDouble);
G3_SERIALIZABLE_CODE(G3VectorString);
G3_SERIALIZABLE_CODE(G3VectorVectorString);
G3_SERIALIZABLE_CODE(G3VectorFrameObject);
G3_SERIALIZABLE_CODE(G3VectorUnsignedChar);
G3_SERIALIZABLE_CODE(G3VectorTime);

/* Special load/save for int64_t. */

static
int bit_count(std::vector<int64_t> const &d) {
	// Returns the smallest number N such that all ints in the
	// vector could be safely expressed as intN_t.  Assumes two's
	// complement integers.  Return value will be between 1 and
	// 64.
	uint64_t mask = 0;
	for (auto c: d) {
		if (c < 0)
			mask |= ~c;
		else
			mask |= c;
	}
	for (int i=1; i<64; i++) {
		if (mask == 0)
			return i;
		mask >>= 1;
	}
	return 64;
}

template <class A, typename FROM_TYPE, typename TO_TYPE>
void load_as(A &ar, std::vector<TO_TYPE> &dest) {
	std::vector<FROM_TYPE> temp;
	ar & cereal::make_nvp("vector", temp);
	dest.resize(temp.size());
	std::copy(temp.begin(), temp.end(), dest.begin());
}

template <class A, typename FROM_TYPE, typename TO_TYPE>
void save_as(A &ar, const std::vector<FROM_TYPE> &src) {
	std::vector<TO_TYPE> temp(src.begin(), src.end());
	ar & cereal::make_nvp("vector", temp);
}

template <>
template <class A>
void G3Vector<int64_t>::load(A &ar, const unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
			      cereal::base_class<G3FrameObject>(this));
	int store_bits = 32;
	if (v >= 2)
		ar & cereal::make_nvp("store_bits", store_bits);

	switch(store_bits) {
	case 64:
		ar & cereal::make_nvp("vector",
				      cereal::base_class<std::vector<int64_t> >(this));
		break;
	case 32:
		load_as<A, int32_t, int64_t>(ar, *this);
		break;
	case 16:
		load_as<A, int16_t, int64_t>(ar, *this);
		break;
	case 8:
		load_as<A, int8_t, int64_t>(ar, *this);
		break;
	}
}

template <>
template <class A>
void G3Vector<int64_t>::save(A &ar, const unsigned v) const
{
	// v == 2
	ar & cereal::make_nvp("G3FrameObject",
			      cereal::base_class<G3FrameObject>(this));
	// Count the interesting bits, and convert to nearest power of 2.
	int sig_bits = bit_count(*this);
	int store_bits = 8;
	while (store_bits < sig_bits)
		store_bits *= 2;
	ar & cereal::make_nvp("store_bits", store_bits);
	switch(store_bits) {
	case 8:
		save_as<A, int64_t, int8_t>(ar, *this);
		break;
	case 16:
		save_as<A, int64_t, int16_t>(ar, *this);
		break;
	case 32:
		save_as<A, int64_t, int32_t>(ar, *this);
		break;
	default:
		ar & cereal::make_nvp("vector",
				      cereal::base_class<std::vector<int64_t> >(this));
	}		
}

template <>
G3VectorBoolPtr
container_from_object(boost::python::object v)
{
	return numpy_container_from_object<G3VectorBool>(v);
}

template <>
G3VectorDoublePtr
container_from_object(boost::python::object v)
{
	return numpy_container_from_object<G3VectorDouble>(v);
}

template <>
G3VectorIntPtr
container_from_object(boost::python::object v)
{
	return numpy_container_from_object<G3VectorInt>(v);
}

template <>
G3VectorComplexDoublePtr
container_from_object(boost::python::object v)
{
	return complex_numpy_container_from_object<G3VectorComplexDouble>(v);
}

template <>
G3VectorTimePtr
container_from_object(boost::python::object v)
{
        return numpy_container_from_object<G3VectorTime>(v);
}

// NB: std::vector<bool> is incompatible with numpy in terms of memory layout
// (stores bits instead of bytes), so no fast path for numpy is provided for
// G3VectorBool. There are ways to make it work read-only if we need it.

static int
G3VectorDouble_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	return pyvector_getbuffer<G3VectorDouble::value_type>(obj, view, flags,
	    "d");
}

static int
G3VectorInt_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	return pyvector_getbuffer<G3VectorInt::value_type>(obj, view, flags,
	    "q");
}

static int
G3VectorComplexDouble_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	return pyvector_getbuffer<G3VectorComplexDouble::value_type>(obj,
	    view, flags, "Zd");
}

static int
G3VectorTime_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	int err;
	G3Time potemkin[2];
	static Py_ssize_t strides = (uintptr_t)&potemkin[1] -
	    (uintptr_t)&potemkin[0];

	err = pyvector_getbuffer<G3VectorTime::value_type>(obj,
	    view, flags, "q");
	if (err != 0)
		return err;

	// pyvector_getbuffer() has set things up so that the elements in
	// the buffer point to the beginning of each G3Time and item size
	// is the size of a G3Time. We need this to be offset to the actual
	// time stamps, maintain the strides, and adjust the item size.
	// Note that offsetof() doesn't work here until C++17, hence the
	// silliness with uintptr_t.
	view->buf = (char *)view->buf + (uintptr_t)&potemkin[0].time -
	    (uintptr_t)&potemkin[0];
	view->itemsize = sizeof(G3TimeStamp);
	view->strides = &strides;
	return 0;
};

static PyBufferProcs vecdouble_bufferprocs;
static PyBufferProcs veccomplexdouble_bufferprocs;
static PyBufferProcs vecint_bufferprocs;
static PyBufferProcs vectime_bufferprocs;

PYBINDINGS("core") {
	boost::python::object vecdouble = register_g3vector<double>(
	    "G3VectorDouble", "Array of floats. Treat as a serializable "
	    "version of numpy.array(dtype=float64). Can be efficiently cast "
	    "to and from numpy arrays.");
	// Add buffer protocol interface
	PyTypeObject *vdclass = (PyTypeObject *)vecdouble.ptr();
	vecdouble_bufferprocs.bf_getbuffer = G3VectorDouble_getbuffer;
	vdclass->tp_as_buffer = &vecdouble_bufferprocs;
#if PY_MAJOR_VERSION < 3
	vdclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif
	boost::python::object veccomplexdouble =
	   register_g3vector<std::complex<double> >(
	    "G3VectorComplexDouble", "Array of complex floats. Treat as a serializable "
	    "version of numpy.array(dtype=complex128). Can be efficiently cast "
	    "to and from numpy arrays.");
	// Add buffer protocol interface
	PyTypeObject *vcclass = (PyTypeObject *)veccomplexdouble.ptr();
	veccomplexdouble_bufferprocs.bf_getbuffer =
	    G3VectorComplexDouble_getbuffer;
	vcclass->tp_as_buffer = &veccomplexdouble_bufferprocs;
#if PY_MAJOR_VERSION < 3
	vcclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif

	boost::python::object vecint = register_g3vector<int64_t>("G3VectorInt",
	    "Array of integers. Treat as a serializable version of "
	    "numpy.array(dtype=int64). Can be efficiently cast to and from "
	    "numpy arrays.");
	// Add buffer protocol interface
	PyTypeObject *viclass = (PyTypeObject *)vecint.ptr();
	vecint_bufferprocs.bf_getbuffer = G3VectorInt_getbuffer,
	viclass->tp_as_buffer = &vecint_bufferprocs;
#if PY_MAJOR_VERSION < 3
	viclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif

	register_g3vector<bool>("G3VectorBool", "List of booleans.");
	register_g3vector<std::string>("G3VectorString", "List of strings.");
	register_vector_of<G3VectorString>("VectorG3VectorString");
	register_g3vector<G3VectorString>("G3VectorVectorString", "List of "
	    "lists of strings.");
	register_g3vector<G3FrameObjectPtr>("G3VectorFrameObject", "List of "
	    "generic frame objects. Can lead to paradoxes; avoid use of this "
	    "class unless you are sure you need it.");
	register_g3vector<uint8_t>("G3VectorUnsignedChar", "List of 8-bit "
	    "integers");
	boost::python::object vectime =
	    register_g3vector<G3Time>("G3VectorTime", "List of times.");
	// Add buffer protocol interface
	PyTypeObject *vtclass = (PyTypeObject *)vectime.ptr();
	vectime_bufferprocs.bf_getbuffer = G3VectorTime_getbuffer,
	vtclass->tp_as_buffer = &vectime_bufferprocs;
#if PY_MAJOR_VERSION < 3
	vtclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif
}

