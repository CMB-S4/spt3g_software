#include <syslog.h>
#include <cstdio>

#include "core/G3Logging.h"
#include "core/G3SimpleLoggers.h"

/*
 * Construct a logger for passing messages to syslog.
 *
 * Arguments:
 *
 * ident : string
 *     Identifier to prepend to syslog messages.
 * facility : integer
 *     Facility type, e.g. LOG_USER.  See syslog(3) for a full list of accepted
 *     values.
 */
G3SyslogLogger::G3SyslogLogger(const std::string& ident, int facility,
  G3LogLevel level)
  : G3Logger(level), syslog_ident(ident), syslog_facility(facility)
{}

void
G3SyslogLogger::Log(G3LogLevel level, const std::string &unit,
    const std::string &file, int line, const std::string &func,
    const std::string &message)
{
	int priority;
	const char* log_description;

	if(LogLevelForUnit(unit) > level)
		return;

	openlog(syslog_ident.empty() ? NULL : syslog_ident.c_str(),
	    LOG_CONS | LOG_PID | LOG_NDELAY, syslog_facility);

	switch(level) {
	case G3LOG_TRACE:
		log_description = "TRACE";
		priority = LOG_DEBUG;
		break;
	case G3LOG_DEBUG:
		log_description = "DEBUG";
		priority = LOG_DEBUG;
		break;
	case G3LOG_INFO:
		log_description = "INFO";
		priority = LOG_INFO;
		break;
	case G3LOG_NOTICE:
		log_description = "NOTICE";
		priority = LOG_NOTICE;
		break;
	case G3LOG_WARN:
		log_description = "WARN";
		priority = LOG_WARNING;
		break;
	case G3LOG_ERROR:
		log_description = "ERROR";
		priority = LOG_ERR;
		break;
	case G3LOG_FATAL:
		log_description = "FATAL";
		priority = LOG_CRIT;
		break;
	default:
		log_description = "UNKNOWN";
		priority = LOG_DEBUG;
		break;
	}

	syslog(priority, "%s (%s): %s (%s:%d in %s)", log_description,
	    unit.c_str(), message.c_str(), file.c_str(), line, func.c_str());

	closelog();
}

