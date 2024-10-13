#include <pybindings.h>

namespace bp = boost::python;

SPT3G_PYTHON_MODULE(core)
{
	bp::docstring_options docopts(true, true, false);

	// Do everything else
	G3ModuleRegistrator::CallRegistrarsFor("core");
}
