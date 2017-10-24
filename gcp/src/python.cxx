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
	bp::enum_<Experiments::Experiment>("Experiments")
		.value("SPT",   Experiments::SPT)
		.value("BK",    Experiments::BK)
		.value("PB",    Experiments::PB)
	;

	G3ModuleRegistrator::CallRegistrarsFor("gcp");
}

