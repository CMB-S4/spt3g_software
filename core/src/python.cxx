#include <pybindings.h>

#include <boost/python.hpp>
#include <boost/preprocessor.hpp>

#include <G3Constants.h>

namespace bp = boost::python;

// The following implements the headerless module registration code
typedef std::map<std::string, std::vector<void (*)()> > module_reg_t;
static module_reg_t *modregs = NULL;

G3ModuleRegistrator::G3ModuleRegistrator(const char *mod, void (*def)())
{
	if (modregs == NULL)
		modregs = new module_reg_t;
	(*modregs)[mod].push_back(def);
}

void G3ModuleRegistrator::CallRegistrarsFor(const char *mod)
{
	for (auto i = (*modregs)[mod].begin(); i != (*modregs)[mod].end(); i++)
		(*i)();
}

G3PythonContext::G3PythonContext(std::string name, bool hold_gil) :
    name_(name), hold_(false), thread_(nullptr)
{
	if (!Py_IsInitialized())
		return;

	if (hold_gil && !PyGILState_Check()) {
		log_debug("%s: Ensuring GIL acquired", name_.c_str());
		gil_ = PyGILState_Ensure();
		hold_ = true;
	} else if (!hold_gil && PyGILState_Check()) {
		log_debug("%s: Saving Python thread state", name_.c_str());
		thread_ = PyEval_SaveThread();
	}
}

G3PythonContext::~G3PythonContext()
{
	if (hold_) {
		log_debug("%s: Releasing GIL", name_.c_str());
		PyGILState_Release(gil_);
		hold_ = false;
	}

	if (!!thread_) {
		log_debug("%s: Restoring Python thread state", name_.c_str());
		PyEval_RestoreThread(thread_);
		thread_ = nullptr;
	}
}

G3PythonInterpreter::G3PythonInterpreter(bool hold_gil) :
    init_(false)
{
	if (!Py_IsInitialized()) {
		log_debug("Initializing");
		Py_Initialize();
#if (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 7)
		PyEval_InitThreads();
#endif
		init_ = true;
	}

	ctx_ = new G3PythonContext("G3PythonInterpreter", hold_gil);
}

G3PythonInterpreter::~G3PythonInterpreter()
{
	if (!!ctx_) {
		delete ctx_;
		ctx_ = nullptr;
	}

	if (init_) {
		log_debug("Finalizing");
		Py_Finalize();
		init_ = false;
	}
}

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
