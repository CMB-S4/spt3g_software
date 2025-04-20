#include <pybindings.h>

// Create a namespace (importable sub-module) within some parent scope
py::object
export_namespace(py::object scope, std::string name)
{
	std::string modname = py::extract<std::string>(scope.attr("__name__") + "." + name);
	py::object mod(py::handle<>(py::borrowed(PyImport_AddModule(modname.c_str()))));
	mod.attr("__package__") = scope.attr("__name__");
	scope.attr(name.c_str()) = mod;
	return mod;
}

// The following implements the headerless module registration code
typedef std::map<std::string, std::vector<void (*)()> > module_reg_t;
static module_reg_t *modregs = NULL;

G3ModuleRegistrator::G3ModuleRegistrator(const char *mod, void (*def)())
{
	if (modregs == NULL)
		modregs = new module_reg_t;
	log_debug("Adding registrar for module %s", mod);
	(*modregs)[mod].push_back(def);
}

void G3ModuleRegistrator::CallRegistrarsFor(const char *mod)
{
	for (auto i = (*modregs)[mod].begin(); i != (*modregs)[mod].end(); i++) {
		log_debug("Calling registrar for module %s", mod);
		(*i)();
	}
}

G3PythonContext::G3PythonContext(std::string name, bool hold_gil) :
    name_(name), hold_(false), thread_(nullptr)
{
	if (!Py_IsInitialized())
		return;

	if (hold_gil && !PyGILState_Check()) {
		log_debug("%s: Ensuring GIL acquired", name_.c_str());
		gil_ = PyGILState_Ensure();
		hold_ = true;
	} else if (!hold_gil && PyGILState_Check()) {
		log_debug("%s: Saving Python thread state", name_.c_str());
		thread_ = PyEval_SaveThread();
	}
}

G3PythonContext::~G3PythonContext()
{
	if (hold_) {
		log_debug("%s: Releasing GIL", name_.c_str());
		PyGILState_Release(gil_);
		hold_ = false;
	}

	if (!!thread_) {
		log_debug("%s: Restoring Python thread state", name_.c_str());
		PyEval_RestoreThread(thread_);
		thread_ = nullptr;
	}
}

G3PythonInterpreter::G3PythonInterpreter(bool hold_gil) :
    init_(false)
{
	if (!Py_IsInitialized()) {
		log_debug("Initializing");
		Py_Initialize();
#if (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 7)
		PyEval_InitThreads();
#endif
		init_ = true;
	}

	ctx_ = new G3PythonContext("G3PythonInterpreter", hold_gil);
}

G3PythonInterpreter::~G3PythonInterpreter()
{
	if (!!ctx_) {
		delete ctx_;
		ctx_ = nullptr;
	}

	if (init_) {
		log_debug("Finalizing");
		Py_Finalize();
		init_ = false;
	}
}

static void translate_ValueError(py::value_error const &e)
{
	PyErr_SetString(PyExc_ValueError, e.what());
}

static void translate_IndexError(py::index_error const &e)
{
	PyErr_SetString(PyExc_IndexError, e.what());
}

static void translate_TypeError(py::type_error const &e)
{
	PyErr_SetString(PyExc_TypeError, e.what());
}

static void translate_KeyError(py::key_error const &e)
{
	PyErr_SetString(PyExc_KeyError, e.what());
}

static void translate_BufferError(py::buffer_error const &e)
{
	PyErr_SetString(PyExc_BufferError, e.what());
}

PYBINDINGS("core")
{
	py::register_exception_translator<py::value_error>(&translate_ValueError);
	py::register_exception_translator<py::index_error>(&translate_IndexError);
	py::register_exception_translator<py::type_error>(&translate_TypeError);
	py::register_exception_translator<py::key_error>(&translate_KeyError);
	py::register_exception_translator<py::buffer_error>(&translate_BufferError);
}
