#include <pybindings.h>

namespace bp = boost::python;

SPT3G_PYTHON_MODULE(calibration)
{
	// Disable noise in doc strings
	bp::docstring_options docopts(true, true, false);

	// Python bindings for this module
	G3ModuleRegistrator::CallRegistrarsFor("calibration");
}

