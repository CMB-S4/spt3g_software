#include <pybindings.h>
#include <serialization.h>
#include <container_pybindings.h>
#include <G3Vector.h>

G3_SERIALIZABLE_CODE(G3VectorInt);
G3_SERIALIZABLE_CODE(G3VectorDouble);
G3_SERIALIZABLE_CODE(G3VectorString);
G3_SERIALIZABLE_CODE(G3VectorFrameObject);
G3_SERIALIZABLE_CODE(G3VectorUnsignedChar);
G3_SERIALIZABLE_CODE(G3VectorTime);

template <>
G3VectorDoublePtr
container_from_object(boost::python::object v)
{
	G3VectorDoublePtr x(new (G3VectorDouble));
	Py_buffer view;
	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) != -1) {
		if (strcmp(view.format, "d") == 0) {
			x->insert(x->begin(), (double *)view.buf,
			(double *)view.buf + view.len/sizeof(double));
		} else if (strcmp(view.format, "f") == 0) {
			x->resize(view.len/sizeof(float));
			for (size_t i = 0; i < view.len/sizeof(float); i++)
				(*x)[i] = ((float *)view.buf)[i];
		} else if (strcmp(view.format, "i") == 0) {
			x->resize(view.len/sizeof(int));
			for (size_t i = 0; i < view.len/sizeof(int); i++)
				(*x)[i] = ((int *)view.buf)[i];
		} else if (strcmp(view.format, "I") == 0) {
			x->resize(view.len/sizeof(int));
			for (size_t i = 0; i < view.len/sizeof(int); i++)
				(*x)[i] = ((unsigned int *)view.buf)[i];
		} else if (strcmp(view.format, "l") == 0) {
			x->resize(view.len/sizeof(long));
			for (size_t i = 0; i < view.len/sizeof(long); i++)
				(*x)[i] = ((unsigned long *)view.buf)[i];
		} else {
			// We could add more types, but why do that?
			// Let Python do the work for obscure cases
			boost::python::container_utils::extend_container(*x, v);
		}
		PyBuffer_Release(&view);
	} else {
		PyErr_Clear();
		boost::python::container_utils::extend_container(*x, v);
	}

	return x;
}

template <>
G3VectorIntPtr
container_from_object(boost::python::object v)
{
	G3VectorIntPtr x(new (G3VectorInt));
	Py_buffer view;
	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) != -1) {
		if (strcmp(view.format, "i") == 0) {
			x->resize(view.len/sizeof(int));
			for (size_t i = 0; i < view.len/sizeof(int); i++)
				(*x)[i] = ((int *)view.buf)[i];
		} else if (strcmp(view.format, "I") == 0) {
			x->resize(view.len/sizeof(int));
			for (size_t i = 0; i < view.len/sizeof(int); i++)
				(*x)[i] = ((unsigned int *)view.buf)[i];
		} else if (strcmp(view.format, "l") == 0) {
			x->resize(view.len/sizeof(long));
			for (size_t i = 0; i < view.len/sizeof(long); i++)
				(*x)[i] = ((unsigned long *)view.buf)[i];
		} else {
			// We could add more types, but why do that?
			// Let Python do the work for obscure cases
			boost::python::container_utils::extend_container(*x, v);
		}
		PyBuffer_Release(&view);
	} else {
		PyErr_Clear();
		boost::python::container_utils::extend_container(*x, v);
	}

	return x;
}

static int
G3VectorDouble_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	if (view == NULL) {
		PyErr_SetString(PyExc_ValueError, "NULL view");
		return -1;
	}

	view->shape = NULL;

	boost::python::handle<> self(boost::python::borrowed(obj));
	boost::python::object selfobj(self);
	G3VectorDoublePtr ts = boost::python::extract<G3VectorDoublePtr>(selfobj)();
	view->obj = obj;
	view->buf = (void*)&(*ts)[0];
	view->len = ts->size() * sizeof(double);
	view->readonly = 0;
	view->itemsize = sizeof(double);
	if (flags & PyBUF_FORMAT)
		view->format = (char *)"d";
	else
		view->format = NULL;
	view->ndim = 1;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION >= 3)
	// Abuse internal pointer in the absence of smalltable. This is safe
	// on all architectures except MIPS N32.
	view->internal = (void *)ts->size();
	view->shape = (Py_ssize_t *)(&view->internal);
#else
	view->smalltable[0] = ts->size();
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

static int
G3VectorInt_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	if (view == NULL) {
		PyErr_SetString(PyExc_ValueError, "NULL view");
		return -1;
	}

	view->shape = NULL;

	boost::python::handle<> self(boost::python::borrowed(obj));
	boost::python::object selfobj(self);
	G3VectorIntPtr ts = boost::python::extract<G3VectorIntPtr>(selfobj)();
	view->obj = obj;
	view->buf = (void*)&(*ts)[0];
	view->len = ts->size() * sizeof(int32_t);
	view->readonly = 0;
	view->itemsize = sizeof(int32_t);
	if (flags & PyBUF_FORMAT)
		view->format = (char *)"i";
	else
		view->format = NULL;
	view->ndim = 1;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION >= 3)
	// Abuse internal pointer in the absence of smalltable. This is safe
	// on all architectures except MIPS N32.
	view->internal = (void *)ts->size();
	view->shape = (Py_ssize_t *)(&view->internal);
#else
	view->smalltable[0] = ts->size();
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

static PyBufferProcs vecdouble_bufferprocs;
static PyBufferProcs vecint_bufferprocs;

template<>
std::string vec_repr<G3FrameObjectPtr>(boost::python::object self)
{
	using namespace boost::python;
	std::stringstream s;

	s << extract<std::string>(self.attr("__class__").attr("__module__"))()
	    << "." << extract<std::string>(self.attr("__class__").attr("__name__"))()
	    << "([";

	std::vector<G3FrameObjectPtr> &selfobject = extract<std::vector<G3FrameObjectPtr> &>(self)();
	if (selfobject.size() == 1) {
		s << selfobject[0]->Summary();
	} else if (selfobject.size() > 1){
		auto i = selfobject.begin();
		while (i != selfobject.end() - 1) {
			s << (*i)->Summary() << ", ";
			i++;
		}
		s << (*i)->Summary();
	}
	s << "])";

	return s.str();
}

PYBINDINGS("core") {
	boost::python::object vecdouble = register_g3vector<double>(
	    "G3VectorDouble", "Array of floats. Treat as a serializable "
	    "version of numpy.array(dtype=float64). Can be efficiently cast "
	    "to and from numpy arrays.");
	// Add buffer protocol interface
	PyTypeObject *vdclass = (PyTypeObject *)vecdouble.ptr();
	vecdouble_bufferprocs.bf_getbuffer = G3VectorDouble_getbuffer,
	vdclass->tp_as_buffer = &vecdouble_bufferprocs;
#if PY_MAJOR_VERSION < 3
	vdclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
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

	register_g3vector<std::complex<double> >("G3VectorComplexDouble",
	    "Array of complex numbers.");
	register_g3vector<std::string>("G3VectorString", "List of strings.");
	register_g3vector<G3FrameObjectPtr>("G3VectorFrameObject", "List of "
	    "generic frame objects. Can lead to paradoxes; avoid use of this "
	    "class unless you are sure you need it.");
	register_g3vector<uint8_t>("G3VectorUnsignedChar", "List of 8-bit "
	    "integers");
	register_g3vector<G3Time>("G3VectorTime", "List of times.");
}

