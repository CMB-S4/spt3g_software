#include <pybindings.h>
#include <serialization.h>

#include <boost/python.hpp>
#include <boost/preprocessor.hpp>

#include <G3Frame.h>
#include <G3Data.h>
#include <G3Module.h>
#include <G3EventBuilder.h>
#include <G3Pipeline.h>
#include <G3Timestream.h>
#include <G3SimpleLoggers.h>

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

// Use G3Frame::save()/G3Frame::load() to provide a pickle interface for frames
struct g3frame_picklesuite : boost::python::pickle_suite
{
	static boost::python::tuple getstate(boost::python::object obj)
	{
		namespace bp = boost::python;
		std::vector<char> buffer;
		boost::iostreams::filtering_ostream os(
		    boost::iostreams::back_inserter(buffer));
		bp::extract<const G3Frame &>(obj)().save(os);
		os.flush();

		return boost::python::make_tuple(obj.attr("__dict__"),
		    bp::object(bp::handle<>(
		    PyBytes_FromStringAndSize(&buffer[0], buffer.size()))));
	}

	static void setstate(boost::python::object obj,	
	    boost::python::tuple state)
	{
		namespace bp = boost::python;
		Py_buffer view;
		PyObject_GetBuffer(bp::object(state[1]).ptr(), &view,
		    PyBUF_SIMPLE);

		std::vector<char> buf((char *)view.buf, (char *)view.buf + view.len);
		bp::extract<bp::dict>(obj.attr("__dict__"))().update(state[0]);
		bp::extract<G3Frame &>(obj)().load(buf);
		PyBuffer_Release(&view);
	}
};

boost::python::list g3frame_keys(const G3Frame &map)
{
        boost::python::list keys;
	std::vector<std::string> keyvec = map.Keys();

        for (auto i = keyvec.begin(); i != keyvec.end(); i++)
                keys.append(*i);
 
        return keys;
}

static void g3frame_python_put(G3Frame &f, std::string name, bp::object obj)
{
	if (bp::extract<G3FrameObjectPtr>(obj).check())
		f.Put(name, bp::extract<G3FrameObjectPtr>(obj)());
	else if (PyBool_Check(obj.ptr()))
		f.Put(name, boost::make_shared<G3Bool>(bp::extract<bool>(obj)()));
	else if (bp::extract<int64_t>(obj).check())
		f.Put(name, boost::make_shared<G3Int>(bp::extract<int64_t>(obj)()));
	else if (bp::extract<double>(obj).check())
		f.Put(name, boost::make_shared<G3Double>(bp::extract<double>(obj)()));
	else if (bp::extract<std::string>(obj).check())
		f.Put(name, boost::make_shared<G3String>(bp::extract<std::string>(obj)()));
	else {
		PyErr_SetString(PyExc_TypeError, "Object is not a G3FrameObject derivative or a plain-old-data type");
		bp::throw_error_already_set();
	}
}

static bp::object g3frame_python_get(G3Frame &f, std::string name)
{
	// Python doesn't have a concept of const. Add subterfuge.
	G3FrameObjectConstPtr element = f[name];
	if (!element) {
		std::string err = "Key \'" + name + "\' not found";
		PyErr_SetString(PyExc_KeyError, err.c_str());
		bp::throw_error_already_set();
	}

	if (!!boost::dynamic_pointer_cast<const G3Int>(element))
		return bp::object(boost::dynamic_pointer_cast<const G3Int>(element)->value);
	else if (!!boost::dynamic_pointer_cast<const G3Double>(element))
		return bp::object(boost::dynamic_pointer_cast<const G3Double>(element)->value);
	else if (!!boost::dynamic_pointer_cast<const G3String>(element))
		return bp::object(boost::dynamic_pointer_cast<const G3String>(element)->value);
	else if (!!boost::dynamic_pointer_cast<const G3Bool>(element))
		return bp::object(boost::dynamic_pointer_cast<const G3Bool>(element)->value);
	else
		return bp::object(boost::const_pointer_cast<G3FrameObject>(element));
}

static std::string g3frame_str(const G3Frame &f)
{
	std::ostringstream oss;
	oss << f;
	return oss.str();
}

static bp::list g3frame_python_values(G3Frame &f)
{
	bp::list values;
	std::vector<std::string> keyvec = f.Keys();

	for (auto i = keyvec.begin(); i != keyvec.end(); i++)
		values.append(g3frame_python_get(f, *i));

	return values;
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
		} else if (bp::extract<G3FramePtr>(ret).check()) {
			out.push_back(bp::extract<G3FramePtr>(ret)());
		} else if (bp::extract<std::vector<G3FramePtr> >(ret).check()) {
			std::vector<G3FramePtr> outlist =
			    bp::extract<std::vector<G3FramePtr> >(ret)();
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
    (volt)(V)(millivolt)(mV)(microvolt)(uV) \
    \
    (ampere)(amp)(A)(milliamp)(mA)(microamp)(uA)(nanoamp)(nA) \
    \
    (kelvin)(K)(millikelvin)(mK)(microkelvin)(uK)(nanokelvin)(nK) \
    (picokelvin)(pK)(rankine)(R)(snausage) \
    \
    (bar)(b)(millibar)(mb) \
    (byte)(B)(bit)(kilobyte)(KB)(megabyte)(MB)(gigabyte)(GB)

#define UNITS_INTERFACE(r,data,T) \
  static double BOOST_PP_CAT(g3units_return_,T)() { return G3Units::T; }
BOOST_PP_SEQ_FOR_EACH(UNITS_INTERFACE,~,UNITS)
#define G3_UNITS_DEF(r,data,T) \
  .add_static_property(BOOST_PP_STRINGIZE(T), &BOOST_PP_CAT(g3units_return_, T))
struct __XXX_fake_g3units_namespace_XXX {};

// Nonsense boilerplate for POD vector numpy bindings
#define numpy_vector_infrastructure(T, conv) \
template <> \
boost::shared_ptr<std::vector<T> > \
container_from_object(boost::python::object v) \
{ \
	return numpy_container_from_object<std::vector<T> >(v); \
} \
static int \
vector_getbuffer_##T(PyObject *obj, Py_buffer *view, int flags) \
{ \
	return pyvector_getbuffer<T>(obj, view, flags, conv); \
} \
static PyBufferProcs vec_bufferprocs_##T;

#if PY_MAJOR_VERSION < 3
#define numpy_vector_of(T, desc) \
{ \
	boost::python::object cls = register_vector_of<T>(desc); \
	PyTypeObject *vdclass = (PyTypeObject *)cls.ptr(); \
	vec_bufferprocs_##T.bf_getbuffer = vector_getbuffer_##T; \
	vdclass->tp_as_buffer = &vec_bufferprocs_##T; \
	vdclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER; \
}
#else
#define numpy_vector_of(T, desc) \
{ \
	boost::python::object cls = register_vector_of<T>(desc); \
	PyTypeObject *vdclass = (PyTypeObject *)cls.ptr(); \
	vec_bufferprocs_##T.bf_getbuffer = vector_getbuffer_##T; \
	vdclass->tp_as_buffer = &vec_bufferprocs_##T; \
}
#endif

numpy_vector_infrastructure(int32_t, "i")
numpy_vector_infrastructure(double, "d")
numpy_vector_infrastructure(float, "f")

BOOST_PYTHON_MODULE(core)
{
	bp::docstring_options docopts(true, true, false);

	// Units values
	bp::class_<__XXX_fake_g3units_namespace_XXX, boost::noncopyable>(
	  "G3Units",
	    "Units suffixes. If you use these, you don't have to worry about "
	    "units arguments to functions. 1 second is 1*G3Units.s.",
	    bp::no_init)
	      BOOST_PP_SEQ_FOR_EACH(G3_UNITS_DEF,~,UNITS);

	// Some POD types
	register_vector_of<bool>("Bool");
	register_vector_of<int64_t>("Int64");
	register_vector_of<uint64_t>("Uint64");
	numpy_vector_of(int32_t, "Int");
	register_vector_of<uint32_t>("UInt");
	numpy_vector_of(double, "Double");
	register_vector_of<std::complex<double> >("ComplexDouble");
	numpy_vector_of(float, "Float");
	register_vector_of<std::string>("String");
	register_vector_of<unsigned char>("UnsignedChar");
	register_vector_of<G3Time>("G3Time");

	// Internal stuff
	bp::class_<G3FrameObject, G3FrameObjectPtr>("G3FrameObject",
	  "Base class for objects that can be added to a frame. All such "
	  "must inherit from G3FrameObject in C++. Pickle hooks are overridden "
	  "to use the fast internal serialization")
	    .def("Description", &G3FrameObject::Description,
	      "Long-form human-readable description of the object")
	    .def("Summary", &G3FrameObject::Summary,
	      "Short (one-line) description of the object")
	    .def("__str__", &G3FrameObject::Summary)
	    .def_pickle(g3frameobject_picklesuite<G3FrameObject>())
	;
	register_pointer_conversions<G3FrameObject>();

	bp::enum_<G3Frame::FrameType>("G3FrameType")
	    .value("Timepoint",       G3Frame::Timepoint)
	    .value("Housekeeping",    G3Frame::Housekeeping)
	    .value("Observation",     G3Frame::Observation)
	    .value("Scan",            G3Frame::Scan)
	    .value("Map",             G3Frame::Map)
	    .value("InstrumentStatus",G3Frame::InstrumentStatus)
	    .value("PipelineInfo",    G3Frame::PipelineInfo)
	    .value("EndProcessing",   G3Frame::EndProcessing)
	    .value("Calibration",     G3Frame::Calibration)
	    .value("Wiring",          G3Frame::Wiring)
	    .value("GcpSlow",         G3Frame::GcpSlow)
	    .value("None",            G3Frame::None)
	;
	enum_none_converter::from_python<G3Frame::FrameType>();
	register_vector_of<G3Frame::FrameType>("FrameType");

	bp::class_<G3Frame, G3FramePtr>("G3Frame",
	  "Frames are the core datatype of the analysis software. They behave "
	  "like Python dictionaries except that they can only store subclasses "
	  "of core.G3FrameObject and the dictionary keys must be strings. "
	  "Pickling and unpickling them uses internal serialization and is "
	  "extremely fast."
	  "\n\n"
	  "In addition to dictionary-like contents, frames have a type code "
	  "(G3Frame.Type) that designates what kind of data are contained in "
	  "it. These types usually indicates information that changes at "
	  "different rates.")
	    .def(bp::init<G3Frame::FrameType>())
	    .def(bp::init<G3Frame>())
	    .def_readwrite("type", &G3Frame::type, "Type code for frame. "
	      "See general G3Frame docstring.")
	    .def("__setitem__", &g3frame_python_put)
	    .def("__getitem__", &g3frame_python_get)
	    .def("keys", &g3frame_keys, "Returns a list of keys in the frame.")
	    .def("__delitem__", &G3Frame::Delete)
	    .def("values", &g3frame_python_values, "Returns a list of the "
	      "values of the items in the frame.")
	    .def("__contains__", (bool (G3Frame::*)(const std::string &) const)
	       &G3Frame::Has)
	    .def("__str__", &g3frame_str)
	    .def("__len__", &G3Frame::size)
	    .def_pickle(g3frame_picklesuite())
	;
	register_vector_of<G3FramePtr>("Frame");
	register_vector_of<G3FrameObjectPtr>("FrameObject");

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
	      (bp::arg("profile")=false, bp::arg("graph")=false), 
	      "Run pipeline. If profile is True, print execution time "
	      "statistics for each module when complete. If graph is True, "
	      "stores control flow data that can be processed with GraphViz "
	      "once retrieved using GetGraphInfo().")
	    .def("GetGraphInfo", &G3Pipeline::GetGraphInfo,
	      "Get stored control flow information from Run(graph=True)")
	    .def("halt_processing", &G3Pipeline_halt_processing,
	      "Halts all running pipelines after they flush all currently "
	      "in-flight frames. Equivalent to what happens when SIGINT is "
	      "sent to this process. Once set, the first module will not be "
	      "called again.")
	    .staticmethod("halt_processing")
	    .def_readonly("last_frame",
	        &G3Pipeline::last_frame)
	;

	bp::enum_<G3LogLevel>("G3LogLevel")
	    .value("LOG_TRACE",  G3LOG_TRACE)
	    .value("LOG_DEBUG",  G3LOG_DEBUG)
	    .value("LOG_INFO",   G3LOG_INFO)
	    .value("LOG_NOTICE", G3LOG_NOTICE)
	    .value("LOG_WARN",   G3LOG_WARN)
	    .value("LOG_ERROR",  G3LOG_ERROR)
	    .value("LOG_FATAL",  G3LOG_FATAL)
	;

	bp::class_<G3Logger, G3LoggerPtr, boost::noncopyable>("G3Logger",
	  "C++ logging abstract base class", bp::no_init)
	    .add_static_property("global_logger", &GetRootLogger, &SetRootLogger)
	    .def("log", &G3Logger::Log)
	    .def("get_level_for_unit", &G3Logger::LogLevelForUnit)
	    .def("set_level_for_unit", &G3Logger::SetLogLevelForUnit)
	    .def("set_level", &G3Logger::SetLogLevel)
        ;
	register_vector_of<G3LoggerPtr>("G3Logger");

	bp::class_<G3NullLogger, bp::bases<G3Logger>,
	  boost::shared_ptr<G3NullLogger>, boost::noncopyable>(
	  "G3NullLogger", "Logger that does not log. Useful if you don't want log messages");
	bp::class_<G3PrintfLogger, bp::bases<G3Logger>,
	  boost::shared_ptr<G3PrintfLogger>, boost::noncopyable>(
	  "G3PrintfLogger", "Logger that prints error messages to stderr (in color, if stderr is a tty).",
	  bp::init<bp::optional<G3LogLevel> >())
	    .def_readwrite("trim_file_names", &G3PrintfLogger::TrimFileNames)
	    .def_readwrite("timestamps", &G3PrintfLogger::Timestamps)
	;
	bp::class_<G3MultiLogger, bp::bases<G3Logger>,
	  boost::shared_ptr<G3MultiLogger>, boost::noncopyable>(
	  "G3MultiLogger", "Log to multiple loggers at once",
	  bp::init<std::vector<G3LoggerPtr> >());
	bp::class_<G3SyslogLogger, bp::bases<G3Logger>,
	  boost::shared_ptr<G3SyslogLogger>, boost::noncopyable>(
	  "G3SyslogLogger", "Pass log messages to the syslog service. Initialize with "
	  "a string identifier and a logging facility. See syslog(3) for details. Example:\n"
	  "\timport syslog\n\tlogger = core.G3SyslogLogger('myprogram', syslog.LOG_USER)",
	  bp::init<std::string, int, bp::optional<G3LogLevel> >());

	// A few things that need to be in a particular order
	bp::enum_<G3Timestream::TimestreamUnits>("G3TimestreamUnits",
	  "Unit scheme for timestreams and maps. Designates different classes "
	  "of units (power, current, on-sky temperature) rather than choices "
	  "of unit within a class (watts vs. horsepower, or K vs. uK), "
	  "transformations between which are handled by core.G3Units.")
	    .value("None",  G3Timestream::None)
	    .value("Counts",  G3Timestream::Counts)
	    .value("Current",  G3Timestream::Current)
	    .value("Power",  G3Timestream::Power)
	    .value("Resistance",  G3Timestream::Resistance)
	    .value("Tcmb",  G3Timestream::Tcmb)
	;
	enum_none_converter::from_python<G3Timestream::TimestreamUnits>();

	// Do everything else
	G3ModuleRegistrator::CallRegistrarsFor("core");
}

