#ifndef _G3_PYBINDINGS_H
#define _G3_PYBINDINGS_H

#include <G3.h>
#include <G3Frame.h>
#include <G3Logging.h>

#include <boost/python.hpp>
#include <boost/python/suite/indexing/container_utils.hpp>
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <container_conversions.h>
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

#include <std_map_indexing_suite.hpp>

#ifdef __clang__
#pragma GCC diagnostic ignored "-Wself-assign-overloaded"
#endif

template <typename T>
void
register_pointer_conversions()
{
	using boost::python::implicitly_convertible;

	implicitly_convertible<std::shared_ptr<T>, G3FrameObjectPtr>();
	implicitly_convertible<std::shared_ptr<T>, std::shared_ptr<const T> >();
	implicitly_convertible<std::shared_ptr<T>, G3FrameObjectConstPtr>();
}

template <typename T>
static std::shared_ptr<T>
container_from_object(boost::python::object v)
{
	std::shared_ptr<T> x(new T());

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

template <typename T, bool proxy=true>
void
register_map(std::string name, const char *docstring)
{
	using namespace boost::python;
	class_<T, std::shared_ptr<T> >(name.c_str())
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

	SET_LOGGER("G3ModuleRegistrator");
};

#define EXPORT_G3MODULE_AND(mod, T, init, docstring, other_defs)   \
	static void registerfunc##T() { \
		using namespace boost::python; \
		class_<T, bases<G3Module>, std::shared_ptr<T>, \
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
	boost::python::class_<T, boost::python::bases<G3FrameObject>, std::shared_ptr<T> >(#T, docstring, boost::python::initf) \
	    .def(boost::python::init<const T &>()) \
	    .def_pickle(g3frameobject_picklesuite<T>())

#define EXPORT_FRAMEOBJECT_NOINITNAMESPACE(T, initf, docstring) \
	boost::python::class_<T, boost::python::bases<G3FrameObject>, std::shared_ptr<T> >(#T, docstring, initf) \
	    .def(boost::python::init<const T &>()) \
	    .def_pickle(g3frameobject_picklesuite<T>())

// Declare a python module with a name that is a bare token:
//     SPT3G_PYTHON_SUBMODULE(foo, "pkg")
// for a package whose fully qualified name will be pkg.foo
#define SPT3G_PYTHON_SUBMODULE(name, pkg) \
BOOST_PYTHON_MODULE(_lib ## name) { \
	auto mod = boost::python::scope(); \
	mod.attr("__name__") = std::string(pkg) + "." + #name; \
	void (spt3g_init_module_ ## name)(); \
	(spt3g_init_module_ ## name)(); \
} \
void (spt3g_init_module_ ## name)()

// Declare a python module with a name that is a bare token:
//     SPT3G_PYTHON_MODULE(foo)
// for a package whose fully qualified name will be spt3g.foo
#define SPT3G_PYTHON_MODULE(name) \
SPT3G_PYTHON_SUBMODULE(name, "spt3g")

// Create a namespace (importable sub-module) within some parent scope
bp::object export_namespace(bp::object scope, std::string name);

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
