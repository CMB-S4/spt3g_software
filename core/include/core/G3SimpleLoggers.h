#ifndef G3SIMPLELOGGERS_H_INCLUDED
#define G3SIMPLELOGGERS_H_INCLUDED

#include <stdexcept>

#include <G3Logging.h>

class G3NullLogger : public G3Logger {
public:
	virtual void Log(G3LogLevel level, const std::string &unit,
	    const std::string &file, int line, const std::string &func,
	    const std::string &message) {};
};

class G3PrintfLogger : public G3Logger {
public:
	G3PrintfLogger(G3LogLevel default_level = G3DefaultLogLevel);
	virtual void Log(G3LogLevel level, const std::string &unit,
	    const std::string &file, int line, const std::string &func,
	    const std::string &message);

	bool TrimFileNames; // Print only leaf paths (default true)
	bool Timestamps;    // Include timestamps in messages (default false)
private:
	bool tty_;
};

#endif

