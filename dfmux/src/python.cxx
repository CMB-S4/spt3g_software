#include <pybindings.h>

SPT3G_PYTHON_MODULE(dfmux, scope)
{
	// Python bindings dependencies
	py::module_::import("spt3g.core");
}

