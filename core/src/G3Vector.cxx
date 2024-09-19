#include <pybindings.h>
#include <serialization.h>
#include <container_pybindings.h>
#include <G3Vector.h>
#include <complex>
#include "int_storage.h"

// Nonsense boilerplate for POD vector numpy bindings
#define numpy_vector_struct(T, name) \
struct numpy_vector_from_python_##name { \
	numpy_vector_from_python_##name() { \
		boost::python::converter::registry::push_back( \
		    &convertible, &construct, \
		    boost::python::type_id<std::vector<T> >()); \
	} \
	static void *convertible(PyObject* obj_ptr) { \
		Py_buffer view; \
		if (PyObject_GetBuffer(obj_ptr, &view, \
		    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) { \
			PyErr_Clear(); \
			return NULL; \
		} \
		if (view.ndim == 0) { \
			PyBuffer_Release(&view); \
			return NULL; \
		} \
		PyBuffer_Release(&view); \
		return obj_ptr; \
	} \
	static void construct(PyObject* obj_ptr, \
	    boost::python::converter::rvalue_from_python_stage1_data* data) { \
		void* storage = ( \
		    (boost::python::converter::rvalue_from_python_storage<std::vector<T> >*)data)->storage.bytes; \
		new (storage) std::vector<T>; \
		boost::shared_ptr<std::vector<T> > swap_storage = numpy_container_from_object<std::vector<T> >(boost::python::object(boost::python::handle<>(boost::python::borrowed(obj_ptr)))); \
		((std::vector<T> *)(storage))->swap(*swap_storage); \
		data->convertible = storage; \
	} \
};

#define numpy_vector_infrastructure(T, name, conv) \
template <> \
boost::shared_ptr<std::vector<T> > \
container_from_object(boost::python::object v) \
{ \
	return numpy_container_from_object<std::vector<T> >(v); \
} \
static int \
vector_getbuffer_##name(PyObject *obj, Py_buffer *view, int flags) \
{ \
	return pyvector_getbuffer<T>(obj, view, flags, conv); \
} \
static PyBufferProcs vec_bufferprocs_##name; \
numpy_vector_struct(T, name)

#if PY_MAJOR_VERSION < 3
#define numpy_vector_of(T, name, desc) \
{ \
	numpy_vector_from_python_##name(); \
	boost::python::object cls = register_vector_of<T>(desc); \
	PyTypeObject *vdclass = (PyTypeObject *)cls.ptr(); \
	vec_bufferprocs_##name.bf_getbuffer = vector_getbuffer_##name; \
	vdclass->tp_as_buffer = &vec_bufferprocs_##name; \
	vdclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER; \
}
#else
#define numpy_vector_of(T, name, desc) \
{ \
	numpy_vector_from_python_##name(); \
	boost::python::object cls = register_vector_of<T>(desc); \
	PyTypeObject *vdclass = (PyTypeObject *)cls.ptr(); \
	vec_bufferprocs_##name.bf_getbuffer = vector_getbuffer_##name; \
	vdclass->tp_as_buffer = &vec_bufferprocs_##name; \
}
#endif

numpy_vector_infrastructure(int64_t, int64_t, "q")
numpy_vector_infrastructure(uint64_t, uint64_t, "Q")
numpy_vector_infrastructure(int32_t, int32_t, "i")
numpy_vector_infrastructure(uint32_t, uint32_t, "I")
numpy_vector_infrastructure(double, double, "d")
numpy_vector_infrastructure(float, float, "f")

// Apple, for their own insane reasons, defines uint64_t as
// "unsigned long long" even on LP64 systems where longs are
// 64-bit. Because "long long" (not a standard C type!) is not
// actually the same type as "long", even when both are 64-bit
// integers, the uint64_t definition above does not do the right
// thing for size_t on 64-bit Apple systems.
//
// Thanks, Apple. "Think Different!"
#if defined(__APPLE__) && defined(__LP64__)
numpy_vector_struct(size_t, size_t)
numpy_vector_struct(ssize_t, ssize_t)
struct apple_size
{
	static PyObject* convert(const std::vector<size_t> &arg) {
		return boost::python::to_python_value<std::vector<uint64_t> >()(*(std::vector<uint64_t> *)(uintptr_t)(&arg));
	}
};

struct apple_ssize
{
	static PyObject* convert(const std::vector<ssize_t> &arg) {
		return boost::python::to_python_value<std::vector<int64_t> >()(*(std::vector<int64_t> *)(intptr_t)(&arg));
	}
};
#endif

template <> boost::shared_ptr<std::vector<std::complex<float> > >
numpy_container_from_object(boost::python::object v)
{
	return complex_numpy_container_from_object<std::vector<std::complex<float> > >(v);
}
template <> boost::shared_ptr<std::vector<std::complex<double> > >
numpy_container_from_object(boost::python::object v)
{
	return complex_numpy_container_from_object<std::vector<std::complex<double> > >(v);
}

numpy_vector_infrastructure(std::complex<double>, cxdouble, "Zd");
numpy_vector_infrastructure(std::complex<float>, cxfloat, "Zf");

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
		load_as<A, int32_t>(ar, *this);
		break;
	case 16:
		load_as<A, int16_t>(ar, *this);
		break;
	case 8:
		load_as<A, int8_t>(ar, *this);
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
		save_as<A, int8_t>(ar, *this);
		break;
	case 16:
		save_as<A, int16_t>(ar, *this);
		break;
	case 32:
		save_as<A, int32_t>(ar, *this);
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
	view->len = view->itemsize*view->shape[0];
	view->strides = &strides;
	return 0;
};

static PyBufferProcs vecdouble_bufferprocs;
static PyBufferProcs veccomplexdouble_bufferprocs;
static PyBufferProcs vecint_bufferprocs;
static PyBufferProcs vectime_bufferprocs;

PYBINDINGS("core") {
	numpy_vector_of(float, float, "Float");
	numpy_vector_of(double, double, "Double");
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

	numpy_vector_of(std::complex<float>, cxfloat, "ComplexFloat");
	numpy_vector_of(std::complex<double>, cxdouble, "ComplexDouble");
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

	numpy_vector_of(int64_t, int64_t, "Int64");
	numpy_vector_of(uint64_t, uint64_t, "UInt64");
	numpy_vector_of(int32_t, int32_t, "Int");
	numpy_vector_of(uint32_t, uint32_t, "UInt");
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

#if defined(__APPLE__) && defined(__LP64__)
	numpy_vector_from_python_size_t();
	numpy_vector_from_python_ssize_t();
	bp::to_python_converter<std::vector<size_t>, apple_size, false>();
	bp::to_python_converter<std::vector<ssize_t>, apple_ssize, false>();
#endif

	register_vector_of<bool>("Bool");
	register_g3vector<bool>("G3VectorBool", "List of booleans.");

	register_vector_of<std::string>("String");
	register_g3vector<std::string>("G3VectorString", "List of strings.");
	register_vector_of<G3VectorString>("VectorG3VectorString");
	register_g3vector<G3VectorString>("G3VectorVectorString", "List of "
	    "lists of strings.");

	register_g3vector<G3FrameObjectPtr>("G3VectorFrameObject", "List of "
	    "generic frame objects. Can lead to paradoxes; avoid use of this "
	    "class unless you are sure you need it.");

	register_vector_of<unsigned char>("UnsignedChar");
	register_g3vector<uint8_t>("G3VectorUnsignedChar", "List of 8-bit "
	    "integers");

	register_vector_of<G3Time>("G3Time");
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

