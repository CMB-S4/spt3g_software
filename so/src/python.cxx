#define NO_IMPORT_ARRAY

#include <pybindings.h>
#include <so/G3SuperTimestream.h>

// See this header file for discussion of numpy compilation issues.
#include <so/so3g_numpy.h>

namespace bp = boost::python;


SPT3G_PYTHON_MODULE(so)
{
    // Python bindings dependencies
    bp::import("spt3g.core");

    bp::docstring_options docopts(true, true, false);
    G3ModuleRegistrator::CallRegistrarsFor("so");
}

