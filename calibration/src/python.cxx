#include <pybindings.h>

SPT3G_PYTHON_MODULE(calibration, scope)
{
	// Python bindings dependencies
	py::module_::import("spt3g.core");
}

