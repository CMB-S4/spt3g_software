#ifndef _G3_PYBINDINGS_H
#define _G3_PYBINDINGS_H

#include <G3.h>
#include <G3Frame.h>
#include <G3Logging.h>

#include <boost/python.hpp>
#include <exception>

namespace py = boost::python;

template <typename T>
void
register_pointer_conversions()
{
	py::implicitly_convertible<std::shared_ptr<T>, G3FrameObjectPtr>();
	py::implicitly_convertible<std::shared_ptr<T>, std::shared_ptr<const T> >();
	py::implicitly_convertible<std::shared_ptr<T>, G3FrameObjectConstPtr>();
}

// Tool for exporting enum elements called 'None', reserved in Python
// Invoke by doing enum_none_converter::from_python<T>()
struct enum_none_converter {
	template <typename Enum, Enum NoneValue = Enum::None>
	static void from_python()
	{
		py::converter::registry::push_back(
		    &enum_none_converter::convertible,
		    &enum_none_converter::construct<Enum, NoneValue>,
		    py::type_id<Enum>());
	}

	static void* convertible(PyObject* object)
	{
		return (object == Py_None) ? object : NULL;
	}

	template <typename Enum, Enum NoneValue = Enum::None>
	static void construct(
	    PyObject* object,
	    py::converter::rvalue_from_python_stage1_data* data)
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
		py::class_<T, py::bases<G3Module>, std::shared_ptr<T>, \
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

// Declare a python module with a name that is a bare token:
//     SPT3G_PYTHON_SUBMODULE(foo, "pkg")
// for a package whose fully qualified name will be pkg.foo
#define SPT3G_PYTHON_SUBMODULE(name, pkg) \
BOOST_PYTHON_MODULE(_lib ## name) { \
	auto mod = py::scope(); \
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
py::object export_namespace(py::object scope, std::string name);

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

// pybind11 compatibility functions
namespace boost { namespace python {

class gil_scoped_release : public G3PythonContext
{
public:
	gil_scoped_release() : G3PythonContext("release", false) {};
};

class gil_scoped_acquire : public G3PythonContext
{
public:
	gil_scoped_acquire() : G3PythonContext("acquire", true) {};
};

class scoped_interpreter : public G3PythonInterpreter
{
public:
	scoped_interpreter() : G3PythonInterpreter(true) {};
};

template <typename T>
T cast(const py::object &obj)
{
	return py::extract<T>(obj)();
}

template <typename T>
bool isinstance(const py::object &obj)
{
	return py::extract<T>(obj).check();
}

class value_error : std::exception
{
public:
	std::string msg;
	value_error(std::string msg="") : msg(msg) {};
	const char * what() const throw() { return msg.c_str(); }
};

class index_error : std::exception
{
public:
	std::string msg;
	index_error(std::string msg="") : msg(msg) {};
	const char * what() const throw() { return msg.c_str(); }
};

class type_error : std::exception
{
public:
	std::string msg;
	type_error(std::string msg="") : msg(msg) {};
	const char * what() const throw() { return msg.c_str(); }
};

class key_error : std::exception
{
public:
	std::string msg;
	key_error(std::string msg="") : msg(msg) {};
	const char * what() const throw() { return msg.c_str(); }
};

class buffer_error : std::exception
{
public:
	std::string msg;
	buffer_error(std::string msg="") : msg(msg) {};
	const char * what() const throw() { return msg.c_str(); }
};

}}

#endif
