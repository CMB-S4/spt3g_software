#include <pybindings.h>

#include <string>
#include <exception>
#ifdef __FreeBSD__
#include <sys/endian.h>
#endif

// Extract object names
std::string py_modname(const py::object &obj)
{
	return py::extract<std::string>(obj.attr("__class__").attr("__module__"))();
}

std::string py_objname(const py::object &obj)
{
	return py::extract<std::string>(obj.attr("__class__").attr("__name__"))();
}

std::string py_fullname(const py::object &obj)
{
	return py_modname(obj) + "." + py_objname(obj);
}

std::string check_buffer_format(std::string fmt) {
	// Consume endian definition
	const char *format = &fmt[0];
	if (format[0] == '@' || format[0] == '=')
		format++;
#if BYTE_ORDER == LITTLE_ENDIAN
	else if (format[0] == '<')
		format++;
	else if (format[0] == '>' || format[0] == '!')
		throw py::buffer_error("Does not support big-endian numpy arrays");
#else
	else if (format[0] == '<')
		throw py::buffer_error("Does not support little-endian numpy arrays");
	else if (format[0] == '>' || format[0] == '!')
		format++;
#endif

	return std::string(format);
}

// The following implements the headerless module registration code
typedef std::map<std::string, std::deque<module_reg_func_t> > module_reg_t;
static std::unique_ptr<module_reg_t> modregs;

G3ModuleRegistrator::G3ModuleRegistrator(const char *mod, module_reg_func_t reg)
{
	if (!modregs)
		modregs = std::unique_ptr<module_reg_t>(new module_reg_t);
	log_debug("Adding registrar for module %s", mod);
	(*modregs)[mod].push_back(reg);
}

void G3ModuleRegistrator::CallRegistrarsFor(const char *mod, py::module_ &scope)
{
	auto &regs = (*modregs)[mod];

	// Call registered functions only once
	while (!regs.empty()) {
		auto reg = regs.front();
		regs.pop_front();
		reg(scope);
	}
}

py::object py::module_::def_submodule(const std::string &name) {
	std::string modname = py::extract<std::string>(
	    this->attr("__name__") + "." + name);
	py::object mod(py::handle<>(py::borrowed(
	    PyImport_AddModule(modname.c_str()))));
	mod.attr("__package__") = this->attr("__name__");
	this->attr(name.c_str()) = mod;
	return mod;
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

G3PythonInterpreter::G3PythonInterpreter() :
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
}

G3PythonInterpreter::~G3PythonInterpreter()
{
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

PYBINDINGS("core", scope)
{
	py::register_exception_translator<py::value_error>(&translate_ValueError);
	py::register_exception_translator<py::index_error>(&translate_IndexError);
	py::register_exception_translator<py::type_error>(&translate_TypeError);
	py::register_exception_translator<py::key_error>(&translate_KeyError);
	py::register_exception_translator<py::buffer_error>(&translate_BufferError);
}
