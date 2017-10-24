#include <pybindings.h>
#include <G3Frame.h>
#include <gcp/Experiments.h>

#include <boost/python.hpp>
#include <boost/preprocessor.hpp>

// Allow customization based on different GCP implementations
#define EXPERIMENTS (SPT)(PB)(BK)

#define EXPERIMENTS_INTERFACE(r,data,T) \
  static unsigned BOOST_PP_CAT(experiments_return_,T)() { return Experiments::T; }
BOOST_PP_SEQ_FOR_EACH(EXPERIMENTS_INTERFACE,~,EXPERIMENTS)
#define EXPERIMENTS_DEF(r,data,T) \
  .add_static_property(BOOST_PP_STRINGIZE(T), &BOOST_PP_CAT(experiments_return_, T))
struct __XXX_fake_experiments_namespace_XXX {};

namespace bp = boost::python;

BOOST_PYTHON_MODULE(gcp)
{
	bp::import("spt3g.core");
	bp::docstring_options docopts(true, true, false);

	// Supported Experiments
	bp::class_<__XXX_fake_experiments_namespace_XXX, boost::noncopyable>(
	  "Experiments",
	    "Supported CMB Experiments",
	    bp::no_init)
	      BOOST_PP_SEQ_FOR_EACH(EXPERIMENTS_DEF,~,EXPERIMENTS);

	G3ModuleRegistrator::CallRegistrarsFor("gcp");
}

