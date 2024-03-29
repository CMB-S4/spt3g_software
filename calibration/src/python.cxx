#include <pybindings.h>

namespace bp = boost::python;

SPT3G_PYTHON_MODULE(calibration)
{
	// Python bindings dependencies
	bp::import("spt3g.core");

	// Disable noise in doc strings
	bp::docstring_options docopts(true, true, false);

	// Python bindings for this module
	G3ModuleRegistrator::CallRegistrarsFor("calibration");
}

