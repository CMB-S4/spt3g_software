#include <pybindings.h>
#include <G3Frame.h>
#include <gcp/Experiments.h>

#include <boost/python.hpp>
#include <boost/preprocessor.hpp>

namespace bp = boost::python;

BOOST_PYTHON_MODULE(gcp)
{
	bp::import("spt3g.core");
	bp::docstring_options docopts(true, true, false);

	// Supported Experiments
	bp::enum_<Experiment>("Experiment")
		.value("SPT",   Experiment::SPT)
		.value("BK",    Experiment::BK)
		.value("PB",    Experiment::PB)
	;

	G3ModuleRegistrator::CallRegistrarsFor("gcp");
}

