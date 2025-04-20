#include <pybindings.h>

SPT3G_PYTHON_MODULE(core)
{
	py::docstring_options docopts(true, true, false);

	// Do everything else
	G3ModuleRegistrator::CallRegistrarsFor("core");
}
