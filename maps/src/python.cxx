#include <pybindings.h>

SPT3G_PYTHON_MODULE(maps)
{
	// Python bindings dependencies
	py::import("spt3g.core");
	py::import("spt3g.calibration");

	py::docstring_options docopts(true, true, false);
	G3ModuleRegistrator::CallRegistrarsFor("maps");
}

