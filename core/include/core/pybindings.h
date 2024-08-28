#ifndef _G3_PYBINDINGS_H
#define _G3_PYBINDINGS_H

#include <G3.h>
#include <G3Frame.h>
#include <G3Logging.h>

#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/facilities/overload.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/python.hpp>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#pragma GCC diagnostic pop
#include <boost/python/suite/indexing/container_utils.hpp>

#include <container_conversions.h>
#include <std_map_indexing_suite.hpp>

template <typename T>
void
register_pointer_conversions()
{
	using boost::python::implicitly_convertible;

	implicitly_convertible<boost::shared_ptr<T>, G3FrameObjectPtr>();
	implicitly_convertible<boost::shared_ptr<T>, boost::shared_ptr<const T> >();
	implicitly_convertible<boost::shared_ptr<T>, G3FrameObjectConstPtr>();
}

template <typename T>
static boost::shared_ptr<T>
container_from_object(boost::python::object v)
{
	boost::shared_ptr<T> x(new T());

	boost::python::container_utils::extend_container(*x, v);
	return x;
}

template <typename T>
std::string vec_repr(boost::python::object self)
{
	using namespace boost::python;
	std::stringstream s;

	s << extract<std::string>(self.attr("__class__").attr("__module__"))()
	    << "." << extract<std::string>(self.attr("__class__").attr("__name__"))()
	    << "([";

	extract <std::vector<T> &> extself(self);
	if (extself.check()) {
		std::vector<T> &selfobject = extself();

		int ellip_pos = -1; // Position at which to insert "..."
		if (selfobject.size() > 100)
			ellip_pos = 3;

		if (selfobject.size() > 0)
			s << selfobject[0];
		for (size_t i=1; i<selfobject.size(); ++i) {
			if ((int)i == ellip_pos) {
				s << ", ...";
				i = selfobject.size() - ellip_pos - 1;
			} else
				s << ", " << selfobject[i];
		}
	}
	s << "])";

	return s.str();
}

// Register naive python conversions to and from std::vector<T>.
// See below if you also want efficient numpy conversions.
template <typename T>
boost::python::object
register_vector_of(std::string name)
{
	name += "Vector";
	using namespace boost::python;
	object cls = class_<std::vector<T> >(name.c_str())
	    .def("__init__", make_constructor(container_from_object<std::vector<T> >))
	    .def("__repr__", vec_repr<T>)
	    .def(vector_indexing_suite<std::vector<T>, true>())
	;
	scitbx::boost_python::container_conversions::from_python_sequence<std::vector<T>,  scitbx::boost_python::container_conversions::variable_capacity_policy>();
        return cls;
}

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

	boost::python::handle<> self(boost::python::borrowed(obj));
	boost::python::object selfobj(self);
	boost::python::extract<const std::vector<T> &> ext(selfobj);
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
boost::shared_ptr<T>
numpy_container_from_object(boost::python::object v)
{
	// There's a chance this is actually a copy operation, so try that first
	bp::extract<T &> extv(v);
	if (extv.check())
		return boost::make_shared<T>(extv());

	boost::shared_ptr<T> x(new T);
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
	boost::python::container_utils::extend_container(*x, v);
        
        return x;
}

// Some special handling is needed for complex vectors

template <typename T>
boost::shared_ptr<T>
complex_numpy_container_from_object(boost::python::object v)
{
	boost::shared_ptr<T> x(new T);
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
		boost::python::container_utils::extend_container(*x, v);
	}

	return x;
}

template <typename T, bool proxy=true>
void
register_map(std::string name, const char *docstring)
{
	using namespace boost::python;
	class_<T, boost::shared_ptr<T> >(name.c_str())
	    .def(init<const T &>())
	    .def(std_map_indexing_suite<T, proxy>())
	;
}

// Tool for exporting enum elements called 'None', reserved in Python
// Invoke by doing enum_none_converter::from_python<T>()
struct enum_none_converter {
	template <typename Enum, Enum NoneValue = Enum::None>
	static void from_python()
	{
		boost::python::converter::registry::push_back(
		    &enum_none_converter::convertible,
		    &enum_none_converter::construct<Enum, NoneValue>,
		    boost::python::type_id<Enum>());
	}

	static void* convertible(PyObject* object)
	{
		return (object == Py_None) ? object : NULL;
	}

	template <typename Enum, Enum NoneValue = Enum::None>
	static void construct(
	    PyObject* object,
	    boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		data->convertible = new Enum(NoneValue);
	}
};

class G3ModuleRegistrator {
public:
	G3ModuleRegistrator(const char *mod, void (*def)());
	static void CallRegistrarsFor(const char *mod);
};

#define EXPORT_G3MODULE_AND(mod, T, init, docstring, other_defs)   \
	static void registerfunc##T() { \
		using namespace boost::python; \
		class_<T, bases<G3Module>, boost::shared_ptr<T>, \
		  boost::noncopyable>(#T, docstring, init) \
		    .def_readonly("__g3module__", true) \
                other_defs \
		; \
	} \
	static G3ModuleRegistrator register##T(mod, registerfunc##T);

#define EXPORT_G3MODULE(mod, T, init, docstring) \
    EXPORT_G3MODULE_AND(mod, T, init, docstring, )


#define PYBINDINGS(mod) \
	static void ___pybindings_registerfunc(); \
	static G3ModuleRegistrator ___pybindings_register(mod, ___pybindings_registerfunc); \
	static void ___pybindings_registerfunc() 

#define EXPORT_FRAMEOBJECT(T, initf, docstring) \
	boost::python::class_<T, boost::python::bases<G3FrameObject>, boost::shared_ptr<T> >(#T, docstring, boost::python::initf) \
	    .def(boost::python::init<const T &>()) \
	    .def_pickle(g3frameobject_picklesuite<T>())

#define EXPORT_FRAMEOBJECT_NOINITNAMESPACE(T, initf, docstring) \
	boost::python::class_<T, boost::python::bases<G3FrameObject>, boost::shared_ptr<T> >(#T, docstring, initf) \
	    .def(boost::python::init<const T &>()) \
	    .def_pickle(g3frameobject_picklesuite<T>())

// Declare a python module with a name and the name of its enclosing package scope.
// name should be be a bare token, while pkg should be a string literal, e.g.:
//     SPT3G_PYTHON_MODULE_2(foo, "spt3g.bar")
// for a package whose fully qualified name will be spt3g.bar.foo
#define SPT3G_PYTHON_MODULE_2(name, pkg) \
BOOST_PYTHON_MODULE(name) { \
	namespace bp = boost::python; \
	auto mod = bp::scope(); \
	std::string package_prefix = pkg; \
	std::string full_name = package_prefix + "." + bp::extract<std::string>(mod.attr("__name__"))(); \
	mod.attr("__name__") = full_name; \
	mod.attr("__package__") = package_prefix; \
	void BOOST_PP_CAT(spt3g_init_module_, name)(); \
	BOOST_PP_CAT(spt3g_init_module_, name)(); \
	if(PY_MAJOR_VERSION < 3){ \
		Py_INCREF(mod.ptr()); \
		PyDict_SetItemString(PyImport_GetModuleDict(),full_name.c_str(),mod.ptr()); \
	} \
} \
void BOOST_PP_CAT(spt3g_init_module_, name)()

// Declare a python module with the given name, assuming that the enclosing package 
// is the default "spt3g".
#define SPT3G_PYTHON_MODULE_1(name) SPT3G_PYTHON_MODULE_2(name, "spt3g")

// Declare a python module with a name and optionally the name of its enclosing package scope.
// name should be be a bare token, while if provided the enclosing package name should be a
// string literal.
// If the enclosing package name is not specified, it will default to "spt3g".
#define SPT3G_PYTHON_MODULE(...) BOOST_PP_OVERLOAD(SPT3G_PYTHON_MODULE_,__VA_ARGS__)(__VA_ARGS__)

// Python runtime context to simplify acquiring or releasing the GIL as necessary.
// To use, simply construct the context object where necessary, e.g.
//    G3PythonContext ctx("mycontext", false);
// The context destructor will clean up after itself (releasing the GIL if acquired, and
// vice versa).  If hold_gil is true, the context will ensure the GIL is held at construction,
// and released at destruction.  If hold_gil is false, the context will save the current thread
// state and release the GIL at construction, and re-acquire it at destruction.
class G3PythonContext {
public:
	G3PythonContext(std::string name, bool hold_gil=false);
	~G3PythonContext();

private:
	std::string name_;
	bool hold_;
	PyGILState_STATE gil_;
	PyThreadState *thread_;

	SET_LOGGER("G3PythonContext");
};

// Convenience class for initializing and finalizing the Python interpreter.  This class
// will initialize the python interpeter at construction (e.g. at the beginning of a C++
// compiled program), and immediately initialize the appropriate G3PythonContext depending
// on the value of hold_gil.  At destruction, it will exit the python context and finalize
// the interpreter.  The python interpreter should be initialized only once, typically at
// the beginning of the main program.
class G3PythonInterpreter {
public:
	G3PythonInterpreter(bool hold_gil=false);
	~G3PythonInterpreter();

private:
	bool init_;
	G3PythonContext *ctx_;

	SET_LOGGER("G3PythonInterpreter");
};

#endif
