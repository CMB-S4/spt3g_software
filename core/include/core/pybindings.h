#ifndef _G3_PYBINDINGS_H
#define _G3_PYBINDINGS_H

#include <G3.h>
#include <G3Frame.h>

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
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

	std::vector<T> &selfobject = extract<std::vector<T> &>(self)();
	if (selfobject.size() == 1) {
		s << selfobject[0];
	} else if (selfobject.size() > 1){
		auto i = selfobject.begin();
		while (i != selfobject.end() - 1) {
			s << *i << ", ";
			i++;
		}
		s << *i;
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
	const std::vector<T> &vec =
	    boost::python::extract<const std::vector<T> &>(selfobj)();
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
	boost::shared_ptr<T> x(new T);
	Py_buffer view;
	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) != -1) {
		if (strcmp(view.format, "d") == 0) {
			x->resize(view.len/sizeof(double));
			for (size_t i = 0; i < view.len/sizeof(double); i++)
				(*x)[i] = ((double *)view.buf)[i];
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
	template <typename Enum>
	static void from_python()
	{
		boost::python::converter::registry::push_back(
		    &enum_none_converter::convertible,
		    &enum_none_converter::construct<Enum>,
		    boost::python::type_id<Enum>());
	}

	static void* convertible(PyObject* object)
	{
		return (object == Py_None) ? object : NULL;
	}

	template <typename Enum>
	static void construct(
	    PyObject* object,
	    boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		data->convertible = new Enum(Enum::None);
	}
};

class G3ModuleRegistrator {
public:
	G3ModuleRegistrator(const char *mod, void (*def)());
	static void CallRegistrarsFor(const char *mod);
};

#define EXPORT_G3MODULE(mod, T, init, docstring) \
	static void registerfunc##T() { \
		using namespace boost::python; \
		class_<T, bases<G3Module>, boost::shared_ptr<T>, \
		  boost::noncopyable>(#T, docstring, init) \
		    .def_readonly("__g3module__", true) \
		; \
	} \
	static G3ModuleRegistrator register##T(mod, registerfunc##T);

#define PYBINDINGS(mod) \
	static void ___pybindings_registerfunc(); \
	static G3ModuleRegistrator ___pybindings_register(mod, ___pybindings_registerfunc); \
	static void ___pybindings_registerfunc() 

#define EXPORT_FRAMEOBJECT(T, initf, docstring) \
	boost::python::class_<T, boost::python::bases<G3FrameObject>, boost::shared_ptr<T> >(#T, docstring, boost::python::initf) \
	    .def(boost::python::init<const T &>()) \
	    .def_pickle(g3frameobject_picklesuite<T>())

#endif

