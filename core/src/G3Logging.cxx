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
	char log_message[messagesize + 1];

	va_start(args, format);
	vsprintf(log_message, format, args);

	return std::string(log_message);
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
	char log_message[messagesize + 1];
	sprintf(log_message, "%s (%s): %s (%s:%d in %s)", log_description,
	    unit.c_str(), message.c_str(), file.c_str(), line, func.c_str());

	BasicLog(log_message);
}

void
g3_clogger(G3LogLevel level, const char *unit, const char *file, int line,
    const char *func, const char *format, ...)
{
	va_list args;
	va_start(args, format);

	int messagesize = vsnprintf(NULL, 0, format, args);
	char log_message[messagesize + 1];

	va_start(args, format);
	vsprintf(log_message, format, args);

	GetRootLogger()->Log(level, unit, file, line, func, log_message);
}

