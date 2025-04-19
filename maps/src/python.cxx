#include <pybindings.h>

SPT3G_PYTHON_MODULE(maps, scope)
{
	// Python bindings dependencies
	py::module_::import("spt3g.core");
	py::module_::import("spt3g.calibration");
}

