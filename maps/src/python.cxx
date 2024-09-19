#include <pybindings.h>

namespace bp = boost::python;

SPT3G_PYTHON_MODULE(maps)
{
	// Python bindings dependencies
	bp::import("spt3g._libcore");
	bp::import("spt3g._libcalibration");

	bp::docstring_options docopts(true, true, false);
	G3ModuleRegistrator::CallRegistrarsFor("maps");
}
