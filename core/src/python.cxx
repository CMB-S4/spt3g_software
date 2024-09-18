#include <pybindings.h>
#include <serialization.h>

#include <boost/python.hpp>
#include <boost/preprocessor.hpp>

#include <G3Module.h>
#include <G3EventBuilder.h>
#include <G3Pipeline.h>
#include <G3TimeStamp.h>
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

static bp::list G3Module_Process(G3Module &mod, G3FramePtr ptr)
{
	std::deque<G3FramePtr> queue;
	bp::list vec;

	mod.Process(ptr, queue);

	for (auto i = queue.begin(); i != queue.end(); i++)
		vec.append(*i);

	return vec;
}

// Python Process() has a different API from the C++ version, since we have
// variable return types. It is passed the frame, but returns its output:
// - If it returns None, the input frame is pushed to the output (usual case)
// - If it returns True, the input frame is pushed to the output
// - If it returns False, the input frame is dropped
// - If it returns a frame, that frame is pushed to the output instead of the
//   input
// - If it returns an iterable of frames, that iterable is pushed instead of
//   the input

class G3ModuleWrap : public G3Module, public boost::python::wrapper<G3Module>
{
public:
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out) {
		bp::object ret = this->get_override("Process")(frame);
		if (ret.ptr() == Py_None) {
			out.push_back(frame);
			return;
		}

		bp::extract<G3FramePtr> extframe(ret);
		if (extframe.check()) {
			out.push_back(extframe());
			return;
		}

		bp::extract<std::vector<G3FramePtr> > extvec(ret);
		if (extvec.check()) {
			std::vector<G3FramePtr> outlist = extvec();
			for (auto i = outlist.begin(); i != outlist.end(); i++)
				out.push_back(*i);
		} else if (!!ret) {
			out.push_back(frame);
		} else {
			// If module returns false on an EndProcessing frame,
			// just let it run through. Doing this is always a bug,
			// so overriding it is safe.
			if (frame->type == G3Frame::EndProcessing)
				out.push_back(frame);
		}
	}
};

static void
G3Pipeline_halt_processing()
{
	G3Pipeline::halt_processing = true;
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

// Nonsense boilerplate for POD vector numpy bindings
#define numpy_vector_struct(T, name) \
struct numpy_vector_from_python_##name { \
	numpy_vector_from_python_##name() { \
		boost::python::converter::registry::push_back( \
		    &convertible, &construct, \
		    boost::python::type_id<std::vector<T> >()); \
	} \
	static void *convertible(PyObject* obj_ptr) { \
		Py_buffer view; \
		if (PyObject_GetBuffer(obj_ptr, &view, \
		    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) { \
			PyErr_Clear(); \
			return NULL; \
		} \
		if (view.ndim == 0) { \
			PyBuffer_Release(&view); \
			return NULL; \
		} \
		PyBuffer_Release(&view); \
		return obj_ptr; \
	} \
	static void construct(PyObject* obj_ptr, \
	    boost::python::converter::rvalue_from_python_stage1_data* data) { \
		void* storage = ( \
		    (boost::python::converter::rvalue_from_python_storage<std::vector<T> >*)data)->storage.bytes; \
		new (storage) std::vector<T>; \
		boost::shared_ptr<std::vector<T> > swap_storage = numpy_container_from_object<std::vector<T> >(boost::python::object(boost::python::handle<>(boost::python::borrowed(obj_ptr)))); \
		((std::vector<T> *)(storage))->swap(*swap_storage); \
		data->convertible = storage; \
	} \
};

#define numpy_vector_infrastructure(T, name, conv) \
template <> \
boost::shared_ptr<std::vector<T> > \
container_from_object(boost::python::object v) \
{ \
	return numpy_container_from_object<std::vector<T> >(v); \
} \
static int \
vector_getbuffer_##name(PyObject *obj, Py_buffer *view, int flags) \
{ \
	return pyvector_getbuffer<T>(obj, view, flags, conv); \
} \
static PyBufferProcs vec_bufferprocs_##name; \
numpy_vector_struct(T, name)

#if PY_MAJOR_VERSION < 3
#define numpy_vector_of(T, name, desc) \
{ \
	numpy_vector_from_python_##name(); \
	boost::python::object cls = register_vector_of<T>(desc); \
	PyTypeObject *vdclass = (PyTypeObject *)cls.ptr(); \
	vec_bufferprocs_##name.bf_getbuffer = vector_getbuffer_##name; \
	vdclass->tp_as_buffer = &vec_bufferprocs_##name; \
	vdclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER; \
}
#else
#define numpy_vector_of(T, name, desc) \
{ \
	numpy_vector_from_python_##name(); \
	boost::python::object cls = register_vector_of<T>(desc); \
	PyTypeObject *vdclass = (PyTypeObject *)cls.ptr(); \
	vec_bufferprocs_##name.bf_getbuffer = vector_getbuffer_##name; \
	vdclass->tp_as_buffer = &vec_bufferprocs_##name; \
}
#endif

numpy_vector_infrastructure(int64_t, int64_t, "q")
numpy_vector_infrastructure(uint64_t, uint64_t, "Q")
numpy_vector_infrastructure(int32_t, int32_t, "i")
numpy_vector_infrastructure(uint32_t, uint32_t, "I")
numpy_vector_infrastructure(double, double, "d")
numpy_vector_infrastructure(float, float, "f")

// Apple, for their own insane reasons, defines uint64_t as
// "unsigned long long" even on LP64 systems where longs are
// 64-bit. Because "long long" (not a standard C type!) is not
// actually the same type as "long", even when both are 64-bit
// integers, the uint64_t definition above does not do the right
// thing for size_t on 64-bit Apple systems.
//
// Thanks, Apple. "Think Different!"
#if defined(__APPLE__) && defined(__LP64__)
numpy_vector_struct(size_t, size_t)
numpy_vector_struct(ssize_t, ssize_t)
struct apple_size
{
	static PyObject* convert(const std::vector<size_t> &arg) {
		return boost::python::to_python_value<std::vector<uint64_t> >()(*(std::vector<uint64_t> *)(uintptr_t)(&arg));
	}
};

struct apple_ssize
{
	static PyObject* convert(const std::vector<ssize_t> &arg) {
		return boost::python::to_python_value<std::vector<int64_t> >()(*(std::vector<int64_t> *)(intptr_t)(&arg));
	}
};
#endif

template <> boost::shared_ptr<std::vector<std::complex<float> > >
numpy_container_from_object(boost::python::object v)
{
	return complex_numpy_container_from_object<std::vector<std::complex<float> > >(v);
}
template <> boost::shared_ptr<std::vector<std::complex<double> > >
numpy_container_from_object(boost::python::object v)
{
	return complex_numpy_container_from_object<std::vector<std::complex<double> > >(v);
}

numpy_vector_infrastructure(std::complex<double>, cxdouble, "Zd");
numpy_vector_infrastructure(std::complex<float>, cxfloat, "Zf");

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

	// Some POD types
	register_vector_of<bool>("Bool");
	numpy_vector_of(int64_t, int64_t, "Int64");
	numpy_vector_of(uint64_t, uint64_t, "UInt64");
	numpy_vector_of(int32_t, int32_t, "Int");
	numpy_vector_of(uint32_t, uint32_t, "UInt");

#if defined(__APPLE__) && defined(__LP64__)
	numpy_vector_from_python_size_t();
	numpy_vector_from_python_ssize_t();
	bp::to_python_converter<std::vector<size_t>, apple_size, false>();
	bp::to_python_converter<std::vector<ssize_t>, apple_ssize, false>();
#endif

	numpy_vector_of(double, double, "Double");
	numpy_vector_of(std::complex<double>, cxdouble, "ComplexDouble");
	numpy_vector_of(float, float, "Float");
	numpy_vector_of(std::complex<float>, cxfloat, "ComplexFloat");
	register_vector_of<std::string>("String");
	register_vector_of<unsigned char>("UnsignedChar");
	register_vector_of<G3Time>("G3Time");

	bp::class_<G3ModuleWrap, boost::shared_ptr<G3ModuleWrap>,
	  boost::noncopyable>("G3Module", "Base class for functors that can be "
	  "added to a G3Pipeline.")
	    .def("__call__", &G3Module_Process)
	    .def("Process", boost::python::pure_virtual(&G3Module_Process))
	;
	bp::implicitly_convertible<boost::shared_ptr<G3ModuleWrap>, G3ModulePtr>();

	bp::class_<G3EventBuilder, bp::bases<G3Module>, G3EventBuilderPtr,
	  boost::noncopyable>("G3EventBuilder", bp::no_init)
	;

	bp::class_<G3Pipeline, boost::shared_ptr<G3Pipeline> >("G3Pipeline",
	  "A collection of core.G3Modules and Python callables. Added "
	  "callables are called sequentially and are passed a frame as their "
	  "only positional argument. If the added callable is a python "
	  "function, it is also passed any keyword arguments given to Add(). "
	  "If it is a class, those keyword arguments are passed to the "
	  "class constructor."
	  "\n\n"
	  "The first module is passed None and returns one or more frames to "
	  "be passed to the next module. Processing will halt if it returns []."
	  "\n\n"
	  "Following modules are passed the frames, one at a time, returned by "
	  "the previous module. The return value from these modules then "
	  "becomes the input queue for the next. Once the last module returns, "
	  "or a module returns [] (or False, which is equivalent), control "
	  "returns to the first module and new data is pushed through the pipe."
	  "\n\n"
	  "Return value semantics for modules:\n"
	  "\t- A single frame: pass frame to next module\n"
	  "\t- An iterable (e.g. a list) of frames: pass frames to next module "
	  "\t  in order\n"
	  "\t- None: pass input frame to next module (implicit if no return)\n"
	  "\t- True: pass input frame to next module\n"
	  "\t- False: discard input frame and return to first module, or end "
	  "\t  processing if returned by first module. Equivalent to [].\n")
	    .def("_Add_", &G3Pipeline::Add, bp::arg("name")="")
	    .def("Run", &G3Pipeline::Run,
	      (bp::arg("profile")=false, bp::arg("graph")=false,
	       bp::arg("signal_halt")=true),
	      "Run pipeline. If profile is True, print execution time "
	      "statistics for each module when complete. If graph is True, "
	      "stores control flow data that can be processed with GraphViz "
	      "once retrieved using GetGraphInfo().  If signal_halt is True "
	      "(default), the pipeline will stop processing new frames when "
	      "SIGINT is sent to this process.  Equivalent to what happens when "
	      "halt_processing() is called.")
	    .def("GetGraphInfo", &G3Pipeline::GetGraphInfo,
	      "Get stored control flow information from Run(graph=True)")
	    .def("halt_processing", &G3Pipeline_halt_processing,
	      "Halts all running pipelines after they flush all currently "
	      "in-flight frames. Once set, the first module will not be "
	      "called again.")
	    .staticmethod("halt_processing")
	    .def_readonly("last_frame",
	        &G3Pipeline::last_frame)
	;

	// Do everything else
	G3ModuleRegistrator::CallRegistrarsFor("core");
}

