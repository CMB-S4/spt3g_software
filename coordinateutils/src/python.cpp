#include <pybindings.h>

namespace bp = boost::python;
BOOST_PYTHON_MODULE(coordinateutils)
{
  bp::import("spt3g.core");
  bp::docstring_options docopts(true, true, false);
  G3ModuleRegistrator::CallRegistrarsFor("coordinateutils");
}
