#include <pybindings.h>

namespace bp = boost::python;

SPT3G_PYTHON_MODULE(maps)
{
	bp::docstring_options docopts(true, true, false);
	G3ModuleRegistrator::CallRegistrarsFor("maps");
}

