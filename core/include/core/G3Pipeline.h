#ifndef _G3_PIPELINE_H
#define _G3_PIPELINE_H

#include <G3Module.h>
#include <vector>
#include <string>
#include <signal.h>

/*
 * Core data processing class. Runs a sequence of modules like a for loop.
 * For example,
 *
 * for frame in frames:
 *   thing1(frame)
 *   thing2(frame)
 *
 * is more or less (see G3Module docs for more complex cases) equivalent to:
 * G3Pipeline pipe()
 *   pipe.Add(thing1)
 *   pipe.Add(thing2)
 * pipe.Run()
 */

class G3Pipeline {
public:
	G3Pipeline();

	// Add a module as the next element of the pipeline. The name argument
	// will be printed as part of profiling output. If unspecified, it is
	// set to the name of the class.
	void Add(G3ModulePtr module, std::string name = "");

	// Run the pipeline to completion. If profile is set to true, will print
	// statistics for each module at completion.
	void Run(bool profile = false, bool graph = false);

	// If there is stored graph information get it
	std::string GetGraphInfo() const {
		return graph_info_;
	}

	// If profile is True in a running pipeline, Run() will record the
	// currently-running module name in a global variable, which can
	// be read using GetCurrentModule(). This is provided on a best-effort
	// basis only and is meaningless unless profile is true.
	static std::string GetCurrentModule();

	static volatile bool halt_processing;
	
	//Stores the last frame of processing if there is an exception
	G3FramePtr last_frame;
private:
	std::vector<std::pair<std::string, G3ModulePtr> > modules_;
	std::string graph_info_;

	static void sigint_catcher(int);
#ifdef SIGINFO
	static void siginfo_catcher(int);
#endif

	SET_LOGGER("G3Pipeline");
};

#endif

