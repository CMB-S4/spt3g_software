#include <pybindings.h>
#include <container_pybindings.h>
#include <G3Logging.h>
#include <G3SimpleLoggers.h>

#include <stdarg.h>

static G3LoggerPtr _global_logger;

G3LoggerPtr
GetRootLogger()
{
	if (!_global_logger)
		_global_logger = G3LoggerPtr(new G3PrintfLogger);

	return _global_logger;
}

void
SetRootLogger(G3LoggerPtr logger)
{

        _global_logger = logger;
}

std::string
G3LoggingStringF(const char *format, ...)
{
	va_list args;
	va_start(args, format);

	int messagesize = vsnprintf(NULL, 0, format, args);
	char *log_message = new char[messagesize + 1];

	va_start(args, format);
	vsnprintf(log_message, messagesize + 1, format, args);

	std::string out(log_message);
	delete [] log_message;

	return out;
}

G3Logger::G3Logger(G3LogLevel default_level) :
    default_log_level_(default_level) {}
G3Logger::~G3Logger() {}

G3LogLevel
G3Logger::LogLevelForUnit(const std::string &unit)
{
	std::map<std::string, G3LogLevel>::const_iterator iter =
	    log_levels_.find(unit);
	if (iter == log_levels_.end())
		return default_log_level_;

	return iter->second;
}

void
G3Logger::SetLogLevelForUnit(const std::string &unit, G3LogLevel level)
{
	log_levels_[unit] = level;
}

void
G3Logger::SetLogLevel(G3LogLevel level)
{
	default_log_level_ = level;
}

G3BasicLogger::G3BasicLogger(G3LogLevel level)
    : G3Logger(level) {}

void
G3BasicLogger::Log(G3LogLevel level, const std::string &unit,
    const std::string &file, int line, const std::string &func,
    const std::string &message)
{
	const char *log_description;

	if (LogLevelForUnit(unit) > level)
		return;

	switch (level) {
	case G3LOG_TRACE:
		log_description = "TRACE";
		break;
	case G3LOG_DEBUG:
		log_description = "DEBUG";
		break;
	case G3LOG_INFO:
		log_description = "INFO";
		break;
        case G3LOG_NOTICE:
                log_description = "NOTICE";
                break;
	case G3LOG_WARN:
		log_description = "WARN";
		break;
	case G3LOG_ERROR:
		log_description = "ERROR";
		break;
	case G3LOG_FATAL:
		log_description = "FATAL";
		break;
	default:
		log_description = "UNKNOWN";
		break;
	}

	int messagesize = snprintf(NULL, 0, "%s (%s): %s (%s:%d in %s)",
	    log_description, unit.c_str(), message.c_str(), file.c_str(), line,
	    func.c_str());
	char *log_message = new char[messagesize + 1];
	snprintf(log_message, messagesize + 1, "%s (%s): %s (%s:%d in %s)", log_description,
	    unit.c_str(), message.c_str(), file.c_str(), line, func.c_str());

	BasicLog(log_message);
	delete [] log_message;
}

void
g3_clogger(G3LogLevel level, const char *unit, const char *file, int line,
    const char *func, const char *format, ...)
{
	va_list args;
	va_start(args, format);

	int messagesize = vsnprintf(NULL, 0, format, args);
	char *log_message = new char[messagesize + 1];

	va_start(args, format);
	vsnprintf(log_message, messagesize + 1, format, args);

	GetRootLogger()->Log(level, unit, file, line, func, log_message);
	delete [] log_message;
}

PYBINDINGS("core", scope) {
	register_enum<G3LogLevel>(scope, "G3LogLevel")
	    .value("LOG_TRACE",  G3LOG_TRACE)
	    .value("LOG_DEBUG",  G3LOG_DEBUG)
	    .value("LOG_INFO",   G3LOG_INFO)
	    .value("LOG_NOTICE", G3LOG_NOTICE)
	    .value("LOG_WARN",   G3LOG_WARN)
	    .value("LOG_ERROR",  G3LOG_ERROR)
	    .value("LOG_FATAL",  G3LOG_FATAL)
	;

	register_class_noncopyable<G3Logger>(scope, "G3Logger", "C++ logging abstract base class")
	    .add_static_property("global_logger", &GetRootLogger, &SetRootLogger)
	    .def("log", &G3Logger::Log)
	    .def("get_level_for_unit", &G3Logger::LogLevelForUnit)
	    .def("set_level_for_unit", &G3Logger::SetLogLevelForUnit)
	    .def("set_level", &G3Logger::SetLogLevel)
        ;
	register_vector_of<G3LoggerPtr>("G3Logger");

	register_class_noncopyable<G3NullLogger, G3Logger>(scope, "G3NullLogger",
	    "Logger that does not log. Useful if you don't want log messages");
	register_class_noncopyable<G3PrintfLogger, G3Logger>(scope, "G3PrintfLogger",
	    "Logger that prints error messages to stderr (in color, if stderr is a tty).")
	    .def(py::init<G3LogLevel>((py::arg("default_level")=G3DefaultLogLevel)))
	    .def_readwrite("trim_file_names", &G3PrintfLogger::TrimFileNames)
	    .def_readwrite("timestamps", &G3PrintfLogger::Timestamps)
	;
	register_class_noncopyable<G3MultiLogger, G3Logger>(scope, "G3MultiLogger",
	    "Log to multiple loggers at once")
	    .def(py::init<std::vector<G3LoggerPtr> >())
	;
	register_class<G3SyslogLogger, G3Logger>(scope, "G3SyslogLogger",
	    "Pass log messages to the syslog service. Initialize with a string identifier "
	    "and a logging facility. See syslog(3) for details. Example:\n"
	    "\timport syslog\n\tlogger = core.G3SyslogLogger('myprogram', syslog.LOG_USER)")
	    .def(py::init<std::string, int, G3LogLevel>((py::arg("ident"), py::arg("facility"),
	        py::arg("default_level")=G3DefaultLogLevel)))
	;
}
