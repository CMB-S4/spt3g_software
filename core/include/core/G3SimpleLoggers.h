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

// Aggregate multiple loggers and log to all of them
class G3MultiLogger : public G3Logger {
public:
	G3MultiLogger(std::vector<G3LoggerPtr> loggers) : loggers_(loggers) {}
	virtual void Log(G3LogLevel level, const std::string &unit,
	    const std::string &file, int line, const std::string &func,
	    const std::string &message) {
		for (auto i = loggers_.begin(); i != loggers_.end(); i++)
			(*i)->Log(level, unit, file, line, func, message);
	}
private:
	std::vector<G3LoggerPtr> loggers_;
};

// Pass log messages to syslog
class G3SyslogLogger : public G3Logger {
public:
	G3SyslogLogger(const std::string &ident, int facility,
	    G3LogLevel default_level=G3DefaultLogLevel);

	virtual void Log(G3LogLevel level, const std::string &unit,
	    const std::string &file, int line, const std::string &func,
	    const std::string &message);
private:
	std::string syslog_ident;
	int syslog_facility;
};

#endif

