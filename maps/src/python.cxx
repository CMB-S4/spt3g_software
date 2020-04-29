#include <pybindings.h>
#include <maps/maputils.h>

namespace bp = boost::python;

BOOST_PYTHON_MODULE(maps)
{
	bp::import("spt3g.core");
	bp::docstring_options docopts(true, true, false);
	maputils_pybindings();
	G3ModuleRegistrator::CallRegistrarsFor("maps");
}

