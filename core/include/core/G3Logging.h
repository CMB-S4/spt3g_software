#ifndef G3LOGGING_H_INCLUDED
#define G3LOGGING_H_INCLUDED

/*
 * Logging framework: use log_info()/log_warn()/log_fatal() in your code
 *  instead of printf(). To have class-specific settable log levels, add
 *  SET_LOGGER("classname") somewhere under private:.
 *
 *  This code is 100% taken from bits of IceTray that I wrote.
 */

#include <G3.h>
#include <signal.h>

typedef enum {
	G3LOG_TRACE,
	G3LOG_DEBUG,
	G3LOG_INFO,
	G3LOG_NOTICE,
	G3LOG_WARN,
	G3LOG_ERROR,
	G3LOG_FATAL
} G3LogLevel;

#if defined(__cplusplus)

#ifdef NDEBUG
const G3LogLevel G3DefaultLogLevel = G3LOG_NOTICE;
#else
const G3LogLevel G3DefaultLogLevel = G3LOG_INFO;
#endif

#include <stdexcept>
#include <map>
#include <sstream>

class G3Logger {
public:
	G3Logger(G3LogLevel default_level = G3DefaultLogLevel);
	virtual ~G3Logger();

	virtual void Log(G3LogLevel level, const std::string &unit,
	    const std::string &file, int line, const std::string &func,
	    const std::string &message) = 0;

	virtual G3LogLevel LogLevelForUnit(const std::string &unit);
	virtual void SetLogLevelForUnit(const std::string &unit,
	    G3LogLevel level);

	virtual void SetLogLevel(G3LogLevel level);
private:
	std::map<std::string, G3LogLevel> log_levels_;
	G3LogLevel default_log_level_;
};

class G3BasicLogger : public G3Logger {
public:
	G3BasicLogger(G3LogLevel default_level = G3DefaultLogLevel);

	virtual void Log(G3LogLevel level, const std::string &unit,
	    const std::string &file, int line, const std::string &func,
	    const std::string &message);
	virtual void BasicLog(const std::string &string) = 0;
};

G3_POINTERS(G3Logger);

// Root logger. If class/namespace has a method with the same name in scope,
// it can have its own logger.
G3LoggerPtr GetRootLogger();
void SetRootLogger(G3LoggerPtr);

std::string G3LoggingStringF(const char *format, ...)
    __attribute__((__format__ (__printf__, 1, 2)));
#define G3_LOGGER(level, id, file, line, func, format, ...) \
    GetRootLogger()->Log(level, id, file, line, func, \
    G3LoggingStringF(format, ##__VA_ARGS__))

#define G3_STREAM_LOGGER(level, id, file, line, func, msg, epilogue) \
    do { std::ostringstream _g3_str_logger_str; _g3_str_logger_str << msg; GetRootLogger()->Log(level, \
    id, file, line, func, _g3_str_logger_str.str()); epilogue } while (0)

extern "C" {
#endif // __cplusplus
void g3_clogger(G3LogLevel level, const char *unit, const char *file,
    int line, const char *func, const char *format, ...)
    __attribute__((__format__ (__printf__, 6, 7)));
#ifdef __cplusplus
}
#else
#define G3_LOGGER g3_clogger
#endif

#define SET_LOGGER(X) \
	static inline const char *__logger_id(void) { return X; }

// Set default logger in global namespace
#ifdef PROJECT
#define G3LOGSTR1(X) #X
#define G3LOGSTR(X) G3LOGSTR1(X)
SET_LOGGER(G3LOGSTR(PROJECT));
#undef G3LOGSTR
#undef G3LOGSTR1
#else
SET_LOGGER("Unknown");
#endif

#define log_custom(level, format, ...) G3_LOGGER(level, \
    __logger_id(), __FILE__, __LINE__, __PRETTY_FUNCTION__, format, \
    ##__VA_ARGS__)
#define log_custom_unit(level, unit, format, ...) G3_LOGGER(level, \
    unit, __FILE__, __LINE__, __PRETTY_FUNCTION__, format, \
    ##__VA_ARGS__)

#ifndef G3_COMPILE_OUT_VERBOSE_LOGGING
#define log_trace(format, ...) G3_LOGGER(G3LOG_TRACE, \
    __logger_id(), __FILE__, __LINE__, __PRETTY_FUNCTION__, format, \
    ##__VA_ARGS__)
#define log_debug(format, ...) G3_LOGGER(G3LOG_DEBUG, \
    __logger_id(), __FILE__, __LINE__, __PRETTY_FUNCTION__, format, \
    ##__VA_ARGS__)
#ifdef __cplusplus
#define log_trace_stream(msg) G3_STREAM_LOGGER(G3LOG_TRACE, \
    __logger_id(), __FILE__, __LINE__, __PRETTY_FUNCTION__, msg, )
#define log_debug_stream(msg) G3_STREAM_LOGGER(G3LOG_DEBUG, \
    __logger_id(), __FILE__, __LINE__, __PRETTY_FUNCTION__, msg, )
#endif
#else
#define log_trace(format, ...)
#define log_debug(format, ...)
#define log_trace_stream(msg)
#define log_debug_stream(msg)
#endif

#define log_info(format, ...) G3_LOGGER(G3LOG_INFO, \
    __logger_id(), __FILE__, __LINE__, __PRETTY_FUNCTION__, format, \
    ##__VA_ARGS__)
#define log_notice(format, ...) G3_LOGGER(G3LOG_NOTICE, \
    __logger_id(), __FILE__, __LINE__, __PRETTY_FUNCTION__, format, \
    ##__VA_ARGS__)
#define log_warn(format, ...) G3_LOGGER(G3LOG_WARN, \
    __logger_id(), __FILE__, __LINE__, __PRETTY_FUNCTION__, format, \
    ##__VA_ARGS__)
#define log_error(format, ...) G3_LOGGER(G3LOG_ERROR, \
    __logger_id(), __FILE__, __LINE__, __PRETTY_FUNCTION__, format, \
    ##__VA_ARGS__)

#ifdef __cplusplus
#define log_fatal(format, ...) G3_LOGGER(G3LOG_FATAL, \
    __logger_id(), __FILE__, __LINE__, __PRETTY_FUNCTION__, format, \
    ##__VA_ARGS__), throw std::runtime_error(G3LoggingStringF(format, \
    ##__VA_ARGS__) + " (in " + __PRETTY_FUNCTION__ + ")")
#define log_info_stream(msg) G3_STREAM_LOGGER(G3LOG_INFO, \
    __logger_id(), __FILE__, __LINE__, __PRETTY_FUNCTION__, msg, )
#define log_notice_stream(msg) G3_STREAM_LOGGER(G3LOG_NOTICE, \
    __logger_id(), __FILE__, __LINE__, __PRETTY_FUNCTION__, msg, )
#define log_warn_stream(msg) G3_STREAM_LOGGER(G3LOG_WARN, \
    __logger_id(), __FILE__, __LINE__, __PRETTY_FUNCTION__, msg, )
#define log_error_stream(msg) G3_STREAM_LOGGER(G3LOG_ERROR, \
    __logger_id(), __FILE__, __LINE__, __PRETTY_FUNCTION__, msg, )
#define log_fatal_stream(msg) G3_STREAM_LOGGER(G3LOG_FATAL, \
    __logger_id(), __FILE__, __LINE__, __PRETTY_FUNCTION__, msg, \
    throw std::runtime_error(_g3_str_logger_str.str() + " (in " + __PRETTY_FUNCTION__ + ")");)
#else
#define log_fatal(format, ...) G3_LOGGER(G3LOG_FATAL, \
    __logger_id(), __FILE__, __LINE__, __PRETTY_FUNCTION__, format, \
    ##__VA_ARGS__), kill(getpid(), SIGABRT)
#endif

#define g3_assert(cond) do{ if(!(cond)) log_fatal("Assertion failed: %s", #cond); } while(0)
#ifdef NDEBUG
#define g3_debug_assert(cond) do{ /*nothing*/ } while(0)
#else
#define g3_debug_assert(cond) do{ if(!(cond)) log_fatal("Assertion failed: %s", #cond); } while(0)
#endif

#endif 

