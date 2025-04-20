#include <pybindings.h>

SPT3G_PYTHON_MODULE(calibration)
{
	// Python bindings dependencies
	py::import("spt3g.core");

	// Disable noise in doc strings
	py::docstring_options docopts(true, true, false);

	// Python bindings for this module
	G3ModuleRegistrator::CallRegistrarsFor("calibration");
}

