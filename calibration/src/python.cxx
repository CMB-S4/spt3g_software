#include <pybindings.h>

SPT3G_PYTHON_MODULE(calibration, scope)
{
	// Python bindings dependencies
	scope.import("spt3g.core");
}

