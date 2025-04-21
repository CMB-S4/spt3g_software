#ifndef _G3_PYBINDINGS_H
#define _G3_PYBINDINGS_H

#include <G3.h>
#include <G3Frame.h>
#include <G3Module.h>
#include <G3Logging.h>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/cast.h>

namespace py = pybind11;

/**
 * pybind11 registration infrastructure
 *
 * Each translation unit that contains structures to expose to python includes a
 * binding function created with the PYBINDINGS macro.
 *
 * Modules are created using SPT3G_PYTHON_MODULE, once in the python module
 * translation unit. In the simplest case, the body of this function calls all of
 * the registered functions created with PYBINDINGS in the various translation
 * units in the linked library.
 *
 * All objects are created with std::shared_ptr holders to allow shared memory
 * access with persistent scope, and the dynamic_attr() option to enable
 * duck-typing in Python.
 */

// Register an enum.  Add enum values using the .value() method of this object.
template <typename T, typename... Args>
auto
register_enum(py::module_ &scope, const std::string &name, Args &&... args)
{
	auto cls = py::enum_<T>(scope, name.c_str(), std::forward<Args>(args)...);

	cls.def_property_readonly("real",
	    [](py::object &obj) { return obj.attr("value"); });
	cls.def_property_readonly_static("names",
	    [](py::object &obj) { return obj.attr("__members__"); });
	cls.def_property_readonly_static("values",
	    [](py::object &obj) {
		py::dict out;
		py::dict mem = obj.attr("__members__");
		for (auto item: mem)
			out[item.second.attr("value")] = item.second;
		return out;
	    });

	std::string rname = scope.attr("__name__").cast<std::string>();

	cls.attr("__str__") = py::cpp_function([](py::object &self) {
		return self.attr("name");
	}, py::name("__str__"), py::is_method(cls));

	cls.attr("__repr__") = py::cpp_function(
	    [rname](py::object &self) {
		return py::str("{}.{}.{}").format(rname,
		    self.attr("__class__").attr("__name__"), self.attr("name"));
	    }, py::name("__repr__"), py::is_method(cls));

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
	cls.def(py::init([](const py::none &) { return NoneValue; }));
	py::implicitly_convertible<py::none, T>();

	return cls;
}

// Register a class, with optional base classes.
// Ensure that base classes are registered first to inherit their bound methods.
template <typename T, typename... Bases, typename... Args>
auto
register_class(py::module_ &scope, const std::string &name, Args&&...args)
{
	return py::class_<T, Bases..., std::shared_ptr<T> >(scope, name.c_str(),
	    py::dynamic_attr(), std::forward<Args>(args)...);
}

// Register a G3Module derived class, for inclusion in a G3Pipeline.
template <typename T, typename... Bases, typename... Args>
auto
register_g3module(py::module_ &scope, const std::string &name, Args&&...args)
{
	auto cls = register_class<T, Bases..., G3Module>(scope, name.c_str(),
	    std::forward<Args>(args)...);

	return cls;
}


// Module registration infrastructure.  Binding functions for
// each module are cached by this class, to be called in order
// when the module is created.

// Registration functions take a single scope object as an input argument.
typedef void (*module_reg_func_t)(py::module_ &);

class PYBIND11_EXPORT G3ModuleRegistrator {
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
PYBIND11_MODULE(_lib ## name, scope) { \
	scope.attr("__name__") = std::string(pkg) + "." + #name; \
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
