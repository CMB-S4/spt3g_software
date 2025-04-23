#include <pybindings.h>

SPT3G_PYTHON_MODULE(gcp, scope)
{
	// Python bindings dependencies
	scope.import("spt3g.core");
	scope.import("spt3g.dfmux");
}

