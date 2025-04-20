#include <pybindings.h>

namespace bp = boost::python;

SPT3G_PYTHON_MODULE(gcp)
{
	// Python bindings dependencies
	bp::import("spt3g.core");
	bp::import("spt3g.dfmux");

	bp::docstring_options docopts(true, true, false);

	G3ModuleRegistrator::CallRegistrarsFor("gcp");
}

