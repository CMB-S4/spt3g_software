#include <G3Logging.h>
#include <G3SimpleLoggers.h>

#include <unistd.h>

G3PrintfLogger::G3PrintfLogger(G3LogLevel level)
    : G3Logger(level), TrimFileNames(true), Timestamps(false)
{
	tty_ = isatty(STDERR_FILENO);
}

void
G3PrintfLogger::Log(G3LogLevel level, const std::string &unit,
    const std::string &file, int line, const std::string &func,
    const std::string &message)
{
	const char *log_description;
	const char *log_prolog = "", *file_prolog = "", *log_epilog = "";

	if (LogLevelForUnit(unit) > level)
		return;

	if (tty_) {
		log_prolog = "\x1b[1m";
		file_prolog = "\x1b[1m";
		log_epilog = "\x1b[0m";
	}

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
		if (tty_)
			log_prolog = "\x1b[1;31m";
		break;
	case G3LOG_FATAL:
		log_description = "FATAL";
		if (tty_)
			log_prolog = "\x1b[1;31m";
		break;
	default:
		log_description = "UNKNOWN";
		break;
	}

	std::string trimmed_filename;
	size_t lastslash = file.rfind('/');
	if (lastslash != std::string::npos && TrimFileNames)
		trimmed_filename = file.substr(lastslash+1);
	else
		trimmed_filename = file;

	char timestamp[255];
	memset(timestamp, 0, sizeof(timestamp));
	if (Timestamps) {
		struct tm tm;
		time_t t = time(NULL);
		localtime_r(&t, &tm);
		strftime(timestamp, sizeof(timestamp),
		    " %d-%b-%Y:%H:%M:%S %Z", &tm);
	}

	int messagesize = snprintf(NULL, 0,
	    "%s%s (%s)%s:%s %s (%s%s:%d%s in %s%s%s)\n",
	    log_prolog, log_description, unit.c_str(), timestamp, log_epilog,
	    message.c_str(), file_prolog, trimmed_filename.c_str(), line,
	    log_epilog, file_prolog, func.c_str(), log_epilog);

	char log_message[messagesize + 1];

	sprintf(log_message,
	    "%s%s (%s)%s:%s %s (%s%s:%d%s in %s%s%s)\n",
	    log_prolog, log_description, unit.c_str(), timestamp, log_epilog,
	    message.c_str(), file_prolog, trimmed_filename.c_str(), line,
	    log_epilog, file_prolog, func.c_str(), log_epilog);

	fputs(log_message, stderr);
}

