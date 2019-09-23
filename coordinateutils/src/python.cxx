#include <pybindings.h>
#include <coordinateutils/flatskyprojection.h>
#include <coordinateutils/maputils.h>

namespace bp = boost::python;

BOOST_PYTHON_MODULE(coordinateutils)
{
	bp::import("spt3g.core");
	bp::docstring_options docopts(true, true, false);
	maputils_pybindings();
	G3ModuleRegistrator::CallRegistrarsFor("coordinateutils");
}

