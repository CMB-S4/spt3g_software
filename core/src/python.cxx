#include <pybindings.h>

#include <boost/python.hpp>
#include <boost/preprocessor.hpp>

#include <G3Constants.h>

namespace bp = boost::python;

#define UNITS (nanosecond)(ns)(microsecond)(us)(millisecond)(ms) \
    (second)(sec)(s)(minute)(min)(hour)(h)(day)                    \
    (nanoseconds)(microseconds)(milliseconds)(seconds)(minutes)\
    (hours)(days)(rel)\
    \
    (Hz)(hz)(MHz)(GHz)(kHz) \
    \
    (rad)(radian)(radians)(deg)(degree)(degrees)(arcmin)(arcsec)(rahour)(rahr) \
    (raminute)(rasecond)\
    \
    (meter)(m)(meters)(centimeter)(cm)(millimeter)(mm)(micron)(nanometer)(nm) \
    (kilometer)(km)(au)(AU)(parsec)(pc)(inch)(in)(foot)(ft)		\
    \
    (watt)(W)(milliwatt)(mW)(microwatt)(uW)(nanowatt)(nW)(picowatt)(pW) \
    (attowatt)(aW)(horsepower)(hp) \
    \
    (sr)(steradian)(steradians)(deg2)(sqdeg)(squaredegree)(squaredegrees) \
    (arcmin2)(sqarcmin)(squarearcmin) \
    \
    (jansky)(Jy)(millijansky)(mJy)(megajansky)(MJy) \
    \
    (volt)(V)(millivolt)(mV)(microvolt)(uV) \
    \
    (ampere)(amp)(A)(milliamp)(mA)(microamp)(uA)(nanoamp)(nA) \
    \
    (kelvin)(K)(millikelvin)(mK)(microkelvin)(uK)(nanokelvin)(nK) \
    (picokelvin)(pK)(rankine)(R)(snausage) \
    \
    (bar)(b)(millibar)(mb)(Pa)(kPa) \
    (byte)(B)(bit)(kilobyte)(KB)(megabyte)(MB)(gigabyte)(GB) \
    (gram)(g)(kilogram)(kg)(milligram)(mg)

#define UNITS_INTERFACE(r,data,T) \
  static double BOOST_PP_CAT(g3units_return_,T)() { return G3Units::T; }
BOOST_PP_SEQ_FOR_EACH(UNITS_INTERFACE,~,UNITS)
#define G3_UNITS_DEF(r,data,T) \
  .add_static_property(BOOST_PP_STRINGIZE(T), &BOOST_PP_CAT(g3units_return_, T))
struct __XXX_fake_g3units_namespace_XXX {};

#define CONSTANTS (c)(h)(hbar)(k)(kb)(G)(g0)(e)

#define CONSTANTS_INTERFACE(r,data,T) \
  static double BOOST_PP_CAT(g3constants_return_,T)() { return G3Constants::T; }
BOOST_PP_SEQ_FOR_EACH(CONSTANTS_INTERFACE,~,CONSTANTS)
#define G3_CONSTANTS_DEF(r,data,T) \
  .add_static_property(BOOST_PP_STRINGIZE(T), &BOOST_PP_CAT(g3constants_return_, T))
struct __XXX_fake_g3constants_namespace_XXX {};

SPT3G_PYTHON_MODULE(core)
{
	bp::docstring_options docopts(true, true, false);

	// Units values
	bp::class_<__XXX_fake_g3units_namespace_XXX, boost::noncopyable>(
	  "G3Units",
	    "Units suffixes. If you use these, you don't have to worry about "
	    "units arguments to functions. 1 second is 1*G3Units.s.",
	    bp::no_init)
	      BOOST_PP_SEQ_FOR_EACH(G3_UNITS_DEF,~,UNITS);

	// Constants values
	bp::class_<__XXX_fake_g3constants_namespace_XXX, boost::noncopyable>(
	  "G3Constants",
	    "Mathematical and physical constants in G3Units.",
	    bp::no_init)
	      BOOST_PP_SEQ_FOR_EACH(G3_CONSTANTS_DEF,~,CONSTANTS);

	// Do everything else
	G3ModuleRegistrator::CallRegistrarsFor("core");
}
