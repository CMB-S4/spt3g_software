#include <pybindings.h>

SPT3G_PYTHON_MODULE(maps, scope)
{
	// Python bindings dependencies
	scope.import("spt3g.core");
	scope.import("spt3g.calibration");
}

