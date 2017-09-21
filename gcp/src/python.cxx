#include <pybindings.h>
#include <G3Frame.h>


#include <boost/python.hpp>

namespace bp = boost::python;

BOOST_PYTHON_MODULE(gcp)
{
	bp::import("spt3g.core");
	bp::docstring_options docopts(true, true, false);

	G3ModuleRegistrator::CallRegistrarsFor("gcp");
}

