#include <pybindings.h>
#include <G3Pipeline.h>
#include <G3Data.h>

#include <sys/time.h>
#include <sys/resource.h>
#include <cxxabi.h>
#include <iostream>
#include <mutex>

G3Pipeline::G3Pipeline() :
  dead_(false)
{
	log_trace("Initializing Pipeline");
}

template <typename T>
static std::string
cxx_demangle(const T& v)
{
	int err;
	const char* name = typeid(v).name();
	char *demangled;

	demangled = abi::__cxa_demangle(name, NULL, NULL, &err);

	std::string demangled_name((err == 0) ? demangled : name);
	free(demangled);

	return demangled_name;
}

void
G3Pipeline::Add(G3ModulePtr module, std::string name)
{
	if (name == "")
		name = cxx_demangle(*module);

	log_trace("Adding module \"%s\"", name.c_str());

	modules_.push_back(
	    std::pair<std::string, G3ModulePtr>(name, module));
}

namespace {
struct G3Pipeline_mod_data {
	G3Pipeline_mod_data(std::pair<std::string, G3ModulePtr> &module) {
		frames = 0;
		mod = module.second;
		name = module.first;
		utime.tv_sec = 0; utime.tv_usec = 0; stime = utime;
		maxrss = 0;
	}

	std::string name;
	G3ModulePtr mod;
	int frames;
	struct timeval utime;
	struct timeval stime;
	long maxrss;
	int graph_id;
};

struct G3Pipeline_proc_data{
	G3Pipeline_proc_data(int mod_id, int frame_id, 
			     G3Frame::FrameType frame_type) :
		mod_id(mod_id),
		frame_id(frame_id),
		frame_type(frame_type) {}
	int mod_id;
	int frame_id;
	G3Frame::FrameType frame_type;
};

// Global accounting variables for debugging (see G3Pipeline::GetCurrentModule)
static std::string pipeline_global_mod;
static std::mutex pipeline_global_mod_lock;

// Recursive inner loop of G3Pipeline::Run()
static size_t
PushFrameThroughQueue(G3FramePtr frame, bool profile, bool graph,
    struct rusage &last_rusage, std::vector<G3Pipeline_mod_data> &mods,
    std::vector<G3Pipeline_mod_data>::iterator next_mod,
    int &graph_frame_id_vals, std::deque<G3Pipeline_proc_data> &graph_proc_data,
    G3FramePtr &last_frame)
{
	std::deque<G3FramePtr> outqueue;

	// Records what processing is happening for graphing
	if (graph && !!frame) {
		int local_frame_id;

		if (frame->Has("_G3GraphingFrameId")) {
			local_frame_id = frame->
			    Get<G3Int>("_G3GraphingFrameId")->value;
		} else {
			frame->Put("_G3GraphingFrameId",
			    G3IntConstPtr(new G3Int(graph_frame_id_vals)));
			local_frame_id = graph_frame_id_vals;
			graph_frame_id_vals++;
		}

		graph_proc_data.push_back(G3Pipeline_proc_data(
		    next_mod->graph_id, local_frame_id, frame->type));
	}

	// Record what module is currently running for debugging purposes
	// if profiling is on.
	if (profile) {
		pipeline_global_mod_lock.lock();
		pipeline_global_mod = next_mod->name;
		pipeline_global_mod_lock.unlock();
	}

	// Actually perform the processing
	try {
		log_trace("Pushing frame through module \"%s\"",
		    next_mod->name.c_str());
		next_mod->mod->Process(frame, outqueue);
	} catch (const std::exception &e) {
		log_warn("Exception in module \"%s\" (%s): %s",
		    next_mod->name.c_str(), typeid(e).name(), e.what());
		last_frame = frame;
		throw;
	} catch (...) {
		log_warn("Exception in module \"%s\"", next_mod->name.c_str());
		last_frame = frame;
		throw;
	}

	// Make sure end processing frames are not eaten, since that can
	// cause strange behavior (unclosed output files and the like).
	if (frame && frame->type == G3Frame::EndProcessing) {
		if (outqueue.size() == 0)
			log_fatal("No output on EndProcessing frame in module "
			    "\"%s\"", next_mod->name.c_str());
		if (outqueue.back()->type != G3Frame::EndProcessing)
			log_fatal("Last queued output frame from module "
			    "\"%s\" on EndProcessing not an EndProcessing "
			    "frame.", next_mod->name.c_str());
	}

	if (profile) {
		struct rusage rusage;
		struct timeval deltat;

#ifdef __APPLE__
		getrusage(RUSAGE_SELF, &rusage);
#else
		getrusage(RUSAGE_THREAD, &rusage);
#endif

		timersub(&rusage.ru_utime, &last_rusage.ru_utime, &deltat);
		timeradd(&next_mod->utime, &deltat, &next_mod->utime);

		timersub(&rusage.ru_stime, &last_rusage.ru_stime, &deltat);
		timeradd(&next_mod->stime, &deltat, &next_mod->stime);

		// If memory usage has increased (allow for lazy deallocation),
		// mark the peak RAM while this module was running.
		if (rusage.ru_maxrss > last_rusage.ru_maxrss + 10*1024 /*10MB*/)
			next_mod->maxrss = rusage.ru_maxrss;

		next_mod->frames++;
		last_rusage = rusage;
	}

	// Bail if this is the last module
	if ((next_mod + 1) == mods.end())
		return outqueue.size();

	// Send each frame on recursively
	for (auto i = outqueue.begin(); i != outqueue.end(); i++)
		PushFrameThroughQueue(*i, profile, graph, last_rusage, mods,
		    next_mod + 1, graph_frame_id_vals, graph_proc_data, 
		    last_frame);

	return outqueue.size();
}
}

std::string
G3Pipeline::GetCurrentModule()
{
	std::string rv;

	pipeline_global_mod_lock.lock();
	rv = pipeline_global_mod;
	pipeline_global_mod_lock.unlock();

	return rv;
}

volatile bool G3Pipeline::halt_processing = false;

void
G3Pipeline::sigint_catcher(int)
{
	log_warn("SIGINT received: halting data processing after "
	    "current frame. Send SIGINT again to abort processing "
	    "immediately, which may result in corrupt output files.");
	G3Pipeline::halt_processing = true;
}

#ifdef SIGINFO
void
G3Pipeline::siginfo_catcher(int)
{
	std::string module = G3Pipeline::GetCurrentModule();

	log_notice("SIGINFO received: currently executing module %s.",
	    module.c_str());
}
#endif

void
G3Pipeline::Run(bool profile, bool graph, bool signal_halt)
{
	struct rusage last_rusage;
	std::vector<G3Pipeline_mod_data> mods;
	struct sigaction sigint_catcher, oldsigint;
#ifdef SIGINFO
	struct sigaction siginfo_catcher, oldsiginfo;
#endif

	// Pipelines can only be run once
	if (dead_)
		log_fatal("Pipeline has already been run and needs to be "
		    "re-initialized to run again.");
	
	// Variables needed for graphing processing chain
	int graph_frame_id_vals = 0;
	std::deque<G3Pipeline_proc_data> graph_proc_data;
	if (modules_.size() == 0)
		log_fatal("No mods present when running pipeline");

	if (profile)
#ifdef __APPLE__
		getrusage(RUSAGE_SELF, &last_rusage);
#else
		getrusage(RUSAGE_THREAD, &last_rusage);
#endif

	if (signal_halt) {
		// Catch SIGINT
		sigint_catcher.sa_handler = &G3Pipeline::sigint_catcher;
		sigint_catcher.sa_flags = SA_RESTART | SA_RESETHAND;
		sigemptyset(&sigint_catcher.sa_mask);
		sigaddset(&sigint_catcher.sa_mask, SIGINT);
		sigaction(SIGINT, &sigint_catcher, &oldsigint);
	}

#ifdef SIGINFO
	if (profile) {
		siginfo_catcher.sa_handler = &G3Pipeline::siginfo_catcher;
		siginfo_catcher.sa_flags = SA_RESTART;
		sigemptyset(&siginfo_catcher.sa_mask);
		sigaddset(&siginfo_catcher.sa_mask, SIGINFO);
		sigaction(SIGINFO, &siginfo_catcher, &oldsiginfo);
	}
#endif

	for (auto i = modules_.begin(); i != modules_.end(); i++)
		mods.push_back(G3Pipeline_mod_data(*i));

	if (graph) {
		for (size_t i = 0; i < mods.size(); i++)
			mods[i].graph_id = i;
	}

	// Start with a NULL frame on the first module, stop when it yields no
	// frames
	while (!G3Pipeline::halt_processing && PushFrameThroughQueue(
	    G3FramePtr(), profile, graph, last_rusage, mods, mods.begin(),
	    graph_frame_id_vals, graph_proc_data, last_frame) != 0)
	{}

	// One last EndProcessing frame, starting at the second module, so
	// later code can do any necessary cleanup.
	if (mods.size() != 1) {
		PushFrameThroughQueue(
		    G3FramePtr(new G3Frame(G3Frame::EndProcessing)), profile,
		    graph, last_rusage, mods, mods.begin() + 1,
		    graph_frame_id_vals, graph_proc_data, 
		    last_frame);
	}

	// Mark us dead -- no coming back now
	dead_ = true;
	for (auto i = mods.begin(); i != mods.end(); i++)
		i->mod.reset(); // Free actual module references

	// Restore old handler
	G3Pipeline::halt_processing = false;
	if (signal_halt)
		sigaction(SIGINT, &oldsigint, &sigint_catcher);
#ifdef SIGINFO
	if (profile)
		sigaction(SIGINFO, &oldsiginfo, &siginfo_catcher);
#endif

	if (profile) {
		struct timeval utime_total = {0, 0};
		struct timeval stime_total = {0, 0};
		long last_maxrss = 0;
		std::string balloonmod = "";

		printf("Pipeline profiling results:\n");
		for (auto i = mods.begin(); i != mods.end(); i++) {
			printf("%s: %ld.%06ld user, %ld.%06ld system, %d frames"
			    " (%.6f s per input frame)\n", i->name.c_str(),
			    long(i->utime.tv_sec), long(i->utime.tv_usec),
			    long(i->stime.tv_sec), long(i->stime.tv_usec),
			    i->frames,
			    (double(i->utime.tv_sec + i->stime.tv_sec) +
			     double(i->utime.tv_usec + i->stime.tv_usec)/1e6) /
			     i->frames);
			timeradd(&utime_total, &i->utime, &utime_total);
			timeradd(&stime_total, &i->stime, &stime_total);

			if (i->maxrss > last_maxrss) {
				balloonmod = i->name;
				last_maxrss = i->maxrss;
			}
		}
		printf("Total: %ld.%06ld user, %ld.%06ld system\n",
		    long(utime_total.tv_sec), long(utime_total.tv_usec),
		    long(stime_total.tv_sec), long(stime_total.tv_usec));
		printf("Peak memory consumption (%.1f MB) in module %s\n",
		    double(last_maxrss) / 1024, balloonmod.c_str());
	}

	if (graph) {
		std::ostringstream os;

		// Formatted to be json parsable: I am a sucker for json
		os << std::endl << "{ \"Module List\": [" << std::endl;
		for (size_t i = 0; i < mods.size(); i++) {
			os << "["<< mods[i].graph_id << ", " 
				  << "\"" << mods[i].name << "\"]";
			if (i < mods.size()-1) os << "," << std::endl;
			else os << std::endl;
		}
		os << "]," <<std::endl << " \"Processing List\": [" << std::endl;
		for (size_t i = 0; i < graph_proc_data.size(); i++) {
			os <<"[" << i << ", " << graph_proc_data[i].mod_id 
				  << ", "<< graph_proc_data[i].frame_id<<", "
				  << "\""<<graph_proc_data[i].frame_type <<"\"]";
			if (i < graph_proc_data.size() -1 ) os<< "," <<std::endl;
			else os<< std::endl;
		}
		os << "] "<< std::endl <<" }" << std::endl;
		graph_info_ = os.str();
	}
}

static py::list G3Module_Process(G3Module &mod, G3FramePtr ptr)
{
	std::deque<G3FramePtr> queue;
	py::list vec;

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

class G3ModuleWrap : public G3Module, public py::wrapper<G3Module>
{
public:
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out) {
		py::object ret = this->get_override("Process")(frame);
		if (ret.ptr() == Py_None) {
			out.push_back(frame);
			return;
		}

		py::extract<G3FramePtr> extframe(ret);
		if (extframe.check()) {
			out.push_back(extframe());
			return;
		}

		py::extract<std::vector<G3FramePtr> > extvec(ret);
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

PYBINDINGS("core", scope) {
	py::class_<G3ModuleWrap, std::shared_ptr<G3ModuleWrap>,
	  boost::noncopyable>("G3Module", "Base class for functors that can be "
	  "added to a G3Pipeline.")
	    .def("__call__", &G3Module_Process)
	    .def("Process", py::pure_virtual(&G3Module_Process))
	;
	py::implicitly_convertible<std::shared_ptr<G3ModuleWrap>, G3ModulePtr>();

	py::class_<G3Pipeline, std::shared_ptr<G3Pipeline> >("G3Pipeline",
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
	    .def("_Add_", &G3Pipeline::Add, py::arg("name")="")
	    .def("Run", &G3Pipeline::Run,
	      (py::arg("profile")=false, py::arg("graph")=false,
	       py::arg("signal_halt")=true),
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
}
