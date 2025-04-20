#include <pybindings.h>
#include <serialization.h>
#include <container_pybindings.h>
#include <G3Vector.h>
#include <complex>
#include "int_storage.h"

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

// Nonsense boilerplate for POD vector numpy bindings
// Code to get a python buffer (for numpy) from one of our classes. "pyformat"
// is the python type code for the type T.
template <typename T>
int pyvector_getbuffer(PyObject *obj, Py_buffer *view, int flags,
    const char *pyformat)
{
	if (view == NULL) {
		PyErr_SetString(PyExc_ValueError, "NULL view");
		return -1;
	}

	view->shape = NULL;

	py::handle<> self(py::borrowed(obj));
	py::object selfobj(self);
	py::extract<const std::vector<T> &> ext(selfobj);
	if (!ext.check()) {
		PyErr_SetString(PyExc_ValueError, "Invalid vector");
		view->obj = NULL;
		return -1;
	}
	const std::vector<T> &vec = ext();
	view->obj = obj;
	view->buf = (void*)&vec[0];
	view->len = vec.size() * sizeof(T);
	view->readonly = 0;
	view->itemsize = sizeof(T);
	if (flags & PyBUF_FORMAT)
		view->format = (char *)pyformat;
	else
		view->format = NULL;
	view->ndim = 1;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION >= 3)
	// Abuse internal pointer in the absence of smalltable. This is safe
	// on all architectures except MIPS N32.
	view->internal = (void *)vec.size();
	view->shape = (Py_ssize_t *)(&view->internal);
#else
	view->smalltable[0] = vec.size();
	view->shape = &view->smalltable[0];
	view->internal = NULL;
#endif
	view->strides = &view->itemsize;
	view->suboffsets = NULL;

	// Try to hold onto our collective hats. This is still very dangerous if
	// the vector's underlying vector is resized.
	Py_INCREF(obj);

	return 0;
}

// For numeric vectors, instantiate this as container_from_object ahead of
// registering the class.
//
// XXX: Should be automatic via SFINAE, but that cleanup is for later.
template <typename T>
std::shared_ptr<T>
numpy_container_from_object(py::object v)
{
	// There's a chance this is actually a copy operation, so try that first
	py::extract<T &> extv(v);
	if (extv.check())
		return std::make_shared<T>(extv());

	std::shared_ptr<T> x(new T);
	size_t nelem;
	Py_buffer view;

	// See if we can get a numpy-ish buffer protocol interface next
	if (PyObject_GetBuffer(v.ptr(), &view, PyBUF_FORMAT | PyBUF_STRIDES)
	    == -1)
		goto slowpython;

	// Only 1D arrays here
	if (view.ndim != 1)
		goto bufferwasbad;

	if (view.shape != NULL)
		nelem = view.shape[0];
	else
		nelem = view.len/view.itemsize;
	x->resize(nelem);

	// Try to get a contiguous buffer first for common data types, since
	// contiguous-copy is likely to be optimizable to memcpy(). In the
	// non-contiguous case below, even if the copy is contiguous, the
	// compiler won't see that. Note that Python's "contiguous" doesn't
	// always mean exactly "contiguous", so we need to check strides too.
	if (PyBuffer_IsContiguous(&view, 'A')) {
		if (strcmp(view.format, "d") == 0) {
			if (view.strides[0] != sizeof(double))
				goto noncontiguous;
			for (size_t i = 0; i < nelem; i++)
				(*x)[i] = ((double *)view.buf)[i];
			goto goodbuffer;
		}
	}

noncontiguous:
	if (strcmp(view.format, "d") == 0) {
		for (size_t i = 0; i < nelem; i++)
			(*x)[i] = *(double *)((char *)view.buf +
			    i*view.strides[0]);
	} else if (strcmp(view.format, "f") == 0) {
		for (size_t i = 0; i < nelem; i++)
			(*x)[i] = *(float *)((char *)view.buf +
			    i*view.strides[0]);
	} else if (strcmp(view.format, "n") == 0) {
		for (size_t i = 0; i < nelem; i++)
			(*x)[i] = *(ssize_t *)((char *)view.buf +
			    i*view.strides[0]);
	} else if (strcmp(view.format, "N") == 0) {
		for (size_t i = 0; i < nelem; i++)
			(*x)[i] = *(size_t *)((char *)view.buf +
			    i*view.strides[0]);
	} else if (strcmp(view.format, "?") == 0) {
		for (size_t i = 0; i < nelem; i++)
			(*x)[i] = *(bool *)((char *)view.buf +
			    i*view.strides[0]);
	} else if (strcmp(view.format, "i") == 0) {
		for (size_t i = 0; i < nelem; i++)
			(*x)[i] = *(int *)((char *)view.buf +
			    i*view.strides[0]);
	} else if (strcmp(view.format, "I") == 0) {
		for (size_t i = 0; i < nelem; i++)
			(*x)[i] = *(unsigned int *)((char *)view.buf +
			    i*view.strides[0]);
	} else if (strcmp(view.format, "l") == 0) {
		for (size_t i = 0; i < nelem; i++)
			(*x)[i] = *(long *)((char *)view.buf +
			    i*view.strides[0]);
	} else if (strcmp(view.format, "L") == 0) {
		for (size_t i = 0; i < nelem; i++)
			(*x)[i] = *(unsigned long *)((char *)view.buf +
			    i*view.strides[0]);
	} else if (strcmp(view.format, "q") == 0) {
		for (size_t i = 0; i < nelem; i++)
			(*x)[i] = *(int64_t *)((char *)view.buf +
			    i*view.strides[0]);
	} else if (strcmp(view.format, "Q") == 0) {
		for (size_t i = 0; i < nelem; i++)
			(*x)[i] = *(uint64_t *)((char *)view.buf +
			    i*view.strides[0]);
	} else {
		// We could add more types, but why do that?
		// Let Python do the work for obscure cases
		goto bufferwasbad;
	}

goodbuffer:
	PyBuffer_Release(&view);
        return x;

bufferwasbad:
	PyBuffer_Release(&view);

slowpython:
	// Failing all of that, try to do element-by-element conversion
	// using Python iteration. If *that* fails, we just give up and
	// go home (TypeError will be set inside extend_container()).
	PyErr_Clear();
	x->resize(0);
	py::container_utils::extend_container(*x, v);

        return x;
}

// Some special handling is needed for complex vectors

template <typename T>
std::shared_ptr<T>
complex_numpy_container_from_object(py::object v)
{
	std::shared_ptr<T> x(new T);
	Py_buffer view;
	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) != -1) {
		if (strcmp(view.format, "Zd") == 0) {
			x->resize(view.len/sizeof(std::complex<double>));
			for (size_t i = 0; i < view.len/sizeof(std::complex<double>); i++)
				(*x)[i] = ((std::complex<double> *)view.buf)[i];
		} else if (strcmp(view.format, "Zf") == 0) {
			x->resize(view.len/sizeof(std::complex<float>));
			for (size_t i = 0; i < view.len/sizeof(std::complex<float>); i++)
				(*x)[i] = ((std::complex<float> *)view.buf)[i];
		} else {
			// Fall back to scalar case otherwise
			auto scalar =
			   numpy_container_from_object<std::vector<double> >(v);
			x->resize(scalar->size());
			for (size_t i = 0; i < scalar->size(); i++)
				(*x)[i] = (*scalar)[i];
		}
		PyBuffer_Release(&view);
	} else {
		PyErr_Clear();
		py::container_utils::extend_container(*x, v);
	}

	return x;
}

#define numpy_vector_struct(T, name)	 \
struct numpy_vector_from_python_##name { \
	numpy_vector_from_python_##name() { \
		py::converter::registry::push_back( \
		    &convertible, &construct, \
		    py::type_id<std::vector<T> >()); \
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
	    py::converter::rvalue_from_python_stage1_data* data) { \
		void* storage = ( \
		    (py::converter::rvalue_from_python_storage<std::vector<T> >*)data)->storage.bytes; \
		new (storage) std::vector<T>; \
		std::shared_ptr<std::vector<T> > swap_storage = numpy_container_from_object<std::vector<T> >(py::object(py::handle<>(py::borrowed(obj_ptr)))); \
		((std::vector<T> *)(storage))->swap(*swap_storage); \
		data->convertible = storage; \
	} \
};

#define numpy_vector_infrastructure(T, name, conv) \
template <> \
std::shared_ptr<std::vector<T> > \
container_from_object(py::object v) \
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
	py::object cls = register_vector_of<T>(scope, desc); \
	PyTypeObject *vdclass = (PyTypeObject *)cls.ptr(); \
	vec_bufferprocs_##name.bf_getbuffer = vector_getbuffer_##name; \
	vdclass->tp_as_buffer = &vec_bufferprocs_##name; \
	vdclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER; \
}
#else
#define numpy_vector_of(T, name, desc) \
{ \
	numpy_vector_from_python_##name(); \
	py::object cls = register_vector_of<T>(scope, desc); \
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
		return py::to_python_value<std::vector<uint64_t> >()(*(std::vector<uint64_t> *)(uintptr_t)(&arg));
	}
};

struct apple_ssize
{
	static PyObject* convert(const std::vector<ssize_t> &arg) {
		return py::to_python_value<std::vector<int64_t> >()(*(std::vector<int64_t> *)(intptr_t)(&arg));
	}
};
#endif

template <> std::shared_ptr<std::vector<std::complex<float> > >
numpy_container_from_object(py::object v)
{
	return complex_numpy_container_from_object<std::vector<std::complex<float> > >(v);
}
template <> std::shared_ptr<std::vector<std::complex<double> > >
numpy_container_from_object(py::object v)
{
	return complex_numpy_container_from_object<std::vector<std::complex<double> > >(v);
}

numpy_vector_infrastructure(std::complex<double>, cxdouble, "Zd");
numpy_vector_infrastructure(std::complex<float>, cxfloat, "Zf");

template <>
G3VectorBoolPtr
container_from_object(py::object v)
{
	return numpy_container_from_object<G3VectorBool>(v);
}

template <>
G3VectorDoublePtr
container_from_object(py::object v)
{
	return numpy_container_from_object<G3VectorDouble>(v);
}

template <>
G3VectorIntPtr
container_from_object(py::object v)
{
	return numpy_container_from_object<G3VectorInt>(v);
}

template <>
G3VectorComplexDoublePtr
container_from_object(py::object v)
{
	return complex_numpy_container_from_object<G3VectorComplexDouble>(v);
}

template <>
G3VectorTimePtr
container_from_object(py::object v)
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

PYBINDINGS("core", scope) {
	numpy_vector_of(float, float, "Float");
	numpy_vector_of(double, double, "Double");
	py::object vecdouble = register_g3vector<G3VectorDouble>(scope,
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
	py::object veccomplexdouble = register_g3vector<G3VectorComplexDouble>(scope,
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
	py::object vecint = register_g3vector<G3VectorInt>(scope, "G3VectorInt",
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
	py::to_python_converter<std::vector<size_t>, apple_size, false>();
	py::to_python_converter<std::vector<ssize_t>, apple_ssize, false>();
#endif

	register_vector_of<bool>(scope, "Bool");
	register_g3vector<G3VectorBool>(scope, "G3VectorBool", "List of booleans.");

	register_vector_of<std::string>(scope, "String");
	register_g3vector<G3VectorString>(scope, "G3VectorString", "List of strings.");
	register_vector_of<G3VectorString>(scope, "G3VectorString");
	register_g3vector<G3VectorVectorString>(scope, "G3VectorVectorString",
	    "List of lists of strings.");

	register_g3vector<G3VectorFrameObject>(scope, "G3VectorFrameObject",
	    "List of generic frame objects. Can lead to paradoxes; avoid use of "
	    "this class unless you are sure you need it.");

	register_vector_of<unsigned char>(scope, "UnsignedChar");
	register_g3vector<G3VectorUnsignedChar>(scope, "G3VectorUnsignedChar",
	    "List of 8-bit integers");

	register_vector_of<G3Time>(scope, "G3Time");
	py::object vectime =
	  register_g3vector<G3VectorTime>(scope, "G3VectorTime", "List of times.");
	// Add buffer protocol interface
	PyTypeObject *vtclass = (PyTypeObject *)vectime.ptr();
	vectime_bufferprocs.bf_getbuffer = G3VectorTime_getbuffer,
	vtclass->tp_as_buffer = &vectime_bufferprocs;
#if PY_MAJOR_VERSION < 3
	vtclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif
}

