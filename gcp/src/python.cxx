#include <pybindings.h>
#include <gcp/Experiments.h>

namespace bp = boost::python;

SPT3G_PYTHON_MODULE(gcp)
{
	// Python bindings dependencies
	bp::import("spt3g.core");
	bp::import("spt3g.dfmux");

	bp::docstring_options docopts(true, true, false);

	// Supported Experiments
	bp::enum_<Experiment>("Experiment")
		.value("SPT",   Experiment::SPT)
		.value("BK",    Experiment::BK)
		.value("PB",    Experiment::PB)
	;

	G3ModuleRegistrator::CallRegistrarsFor("gcp");
}

