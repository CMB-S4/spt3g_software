#include <pybindings.h>
#include <serialization.h>
#include <container_pybindings.h>
#include <G3Vector.h>
#include <complex>

G3_SERIALIZABLE_CODE(G3VectorInt);
G3_SERIALIZABLE_CODE(G3VectorDouble);
G3_SERIALIZABLE_CODE(G3VectorComplexDouble);
G3_SERIALIZABLE_CODE(G3VectorString);
G3_SERIALIZABLE_CODE(G3VectorVectorString);
G3_SERIALIZABLE_CODE(G3VectorFrameObject);
G3_SERIALIZABLE_CODE(G3VectorUnsignedChar);
G3_SERIALIZABLE_CODE(G3VectorTime);

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
	    "i");
}

static int
G3VectorComplexDouble_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	return pyvector_getbuffer<G3VectorComplexDouble::value_type>(obj,
	    view, flags, "Zd");
}

static PyBufferProcs vecdouble_bufferprocs;
static PyBufferProcs veccomplexdouble_bufferprocs;
static PyBufferProcs vecint_bufferprocs;

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

	boost::python::object vecint = register_g3vector<int32_t>("G3VectorInt",
	    "Array of integers. Treat as a serializable version of "
	    "numpy.array(dtype=int32). Can be efficiently cast to and from "
	    "numpy arrays.");
	// Add buffer protocol interface
	PyTypeObject *viclass = (PyTypeObject *)vecint.ptr();
	vecint_bufferprocs.bf_getbuffer = G3VectorInt_getbuffer,
	viclass->tp_as_buffer = &vecint_bufferprocs;
#if PY_MAJOR_VERSION < 3
	viclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif

	register_g3vector<std::string>("G3VectorString", "List of strings.");
	register_vector_of<G3VectorString>("VectorG3VectorString");
	register_g3vector<G3VectorString>("G3VectorVectorString", "List of "
	    "lists of strings.");
	register_g3vector<G3FrameObjectPtr>("G3VectorFrameObject", "List of "
	    "generic frame objects. Can lead to paradoxes; avoid use of this "
	    "class unless you are sure you need it.");
	register_g3vector<uint8_t>("G3VectorUnsignedChar", "List of 8-bit "
	    "integers");
	register_g3vector<G3Time>("G3VectorTime", "List of times.");
}

