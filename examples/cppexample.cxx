#include <core/pybindings.h>
#include <core/G3.h>
#include <core/G3Pipeline.h>
#include <core/G3Reader.h>
#include <core/G3Writer.h>

/*
 * Example of a small C++ program that is the equivalent of the
 * spt3g-dump command.
 *
 * Demonstrates how to use G3Pipeline and write modules in C++.
 */

class Dump : public G3Module
{
public:
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
	{
		std::cout << *frame << std::endl;
		out.push_back(frame);
	}
};

int
main(int argc, const char **argv)
{
	if (argc < 2) {
		std::cerr << "Too few arguments!" << std::endl;
		return 1;
	}

	// Initialize the python interpreter, and release the GIL.
	// Set the argument to true to instead hold the GIL.
	// Comment this out to disable the interpreter.
	G3PythonInterpreter interp(false);

	G3Pipeline pipe;

	pipe.Add(G3ModulePtr(new G3Reader(argv[1])));
	pipe.Add(G3ModulePtr(new Dump));

	if (argc > 2)
		pipe.Add(G3ModulePtr(new G3Writer(argv[2])));

	pipe.Run();
	
	return 0;
}

