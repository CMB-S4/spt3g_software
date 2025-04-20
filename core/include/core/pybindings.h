#ifndef _G3_PYBINDINGS_H
#define _G3_PYBINDINGS_H

#include <G3.h>
#include <G3Frame.h>
#include <G3Module.h>
#include <G3Logging.h>

#include <pybind11_compat.h>

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

// Register an enum.  Add enum values using the .value() method of this object.
template <typename T, typename... Args>
auto
register_enum(py::module_ &scope, const std::string &name, Args &&... args)
{
	(void) scope;

	auto cls = py::enum_<T>(name.c_str(), std::forward<Args>(args)...);

	return cls;
}

// Register an enum with a particular value that corresponds to python `None`.
// For example:
//     register_enum<G3Frame::FrameType, G3Frame::None>(scope, "G3FrameType")
// Creates a G3FrameType python object.  A None value passed to any function
// in python that takes this type will be translated to G3Frame::None in C++.
template <typename T, T NoneValue, typename... Args>
auto
register_enum(py::module_ &scope, const std::string &name, Args &&... args)
{
	auto cls = register_enum<T>(scope, name, std::forward<Args>(args)...);

	// Convert python None to a specific value
	enum_none_converter::from_python<T, NoneValue>();

	return cls;
}

// Register a class, with optional base classes.
// Ensure that base classes are registered first to inherit their bound methods.
template <typename T, typename... Bases, typename... Args>
auto
register_class(py::module_ &scope, const std::string &name, Args&&...args)
{
	(void) scope;

	return py::class_<T, py::bases<Bases...>, std::shared_ptr<T> >(name.c_str(),
	    std::forward<Args>(args)..., py::no_init);
}

// Register a class, with optional base classes.
// Ensure that base classes are registered first to inherit their bound methods.
template <typename T, typename... Bases, typename... Args>
auto
register_class_noncopyable(py::module_ &scope, const std::string &name, Args&&...args)
{
	(void) scope;

	return py::class_<T, py::bases<Bases...>, std::shared_ptr<T>,
	    boost::noncopyable>(name.c_str(), std::forward<Args>(args)..., py::no_init);
}

// Register a G3Module derived class, for inclusion in a G3Pipeline.
template <typename T, typename... Bases, typename... Args>
auto
register_g3module(py::module_ &scope, const std::string &name, Args&&...args)
{
	auto cls = register_class_noncopyable<T, Bases..., G3Module>(scope, name.c_str(),
	    std::forward<Args>(args)...);

	// mark as a module
	try {
		cls.def_readonly("__g3module__", true);
	} catch (py::error_already_set &e) {
		PyErr_Clear();
	}

	py::implicitly_convertible<std::shared_ptr<T>, G3ModulePtr>();

	return cls;
}


// Module registration infrastructure.  Binding functions for
// each module are cached by this class, to be called in order
// when the module is created.

// Registration functions take a single scope object as an input argument.
typedef void (*module_reg_func_t)(py::module_ &);

class G3ModuleRegistrator {
public:
	G3ModuleRegistrator(const char *mod, module_reg_func_t def);
	static void CallRegistrarsFor(const char *mod, py::module_ &scope);

	SET_LOGGER("G3ModuleRegistrator");
};

// Create and register a binding function for the given module,
// to be called by the registrar when the python module is constructed.
// This macro can be used at most once per translation unit.
// mod : string name of the python module, e.g. "core"
// scope : variable name of the pybind11::module_ object to be populated
//    use this name as the first argument to the various registration functions
#define PYBINDINGS(mod, scope) \
	static void ___pybindings_registerfunc(py::module_ &); \
	static G3ModuleRegistrator ___pybindings_register(mod, ___pybindings_registerfunc); \
	static void ___pybindings_registerfunc(py::module_ &scope)

// Declare a python module with a name that is a bare token:
//     SPT3G_PYTHON_SUBMODULE(foo, "pkg", scope)
// for a package whose fully qualified name will be pkg.foo
#define SPT3G_PYTHON_SUBMODULE(name, pkg, scope) \
BOOST_PYTHON_MODULE(_lib ## name) { \
	auto scope = py::module_(); \
	scope.attr("__name__") = std::string(pkg) + "." + #name; \
	py::docstring_options docopts(true, true, false); \
	void (spt3g_init_module_ ## name)(py::module_ &); \
	(spt3g_init_module_ ## name)(scope); \
	G3ModuleRegistrator::CallRegistrarsFor(#name, scope); \
} \
void (spt3g_init_module_ ## name)([[maybe_unused]] py::module_ &scope)

// Declare a python module with a name that is a bare token:
//     SPT3G_PYTHON_MODULE(foo, scope)
// for a package whose fully qualified name will be spt3g.foo
//
// Use this macro once per module, as a function whose body is
// populated with preliminaries, e.g. importing necessary dependencies.
// All registered binding functions are called after this block.
#define SPT3G_PYTHON_MODULE(name, scope) \
SPT3G_PYTHON_SUBMODULE(name, "spt3g", scope)

#endif
