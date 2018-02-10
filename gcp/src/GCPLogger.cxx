#include <G3Logging.h>
#include <pybindings.h>

#include <deque>
#include <string>
#include <mutex>
#include <condition_variable>
#include <thread>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>

class GCPLogger : public G3Logger {
public:
	GCPLogger(int port, G3LogLevel default_level);
	virtual ~GCPLogger();
	virtual void Log(G3LogLevel level, const std::string &unit,
	    const std::string &file, int line, const std::string &func,
	    const std::string &message);
	bool TrimFileNames;

private:
	int fd_;
	std::deque<std::string> log_deque_;
	std::mutex log_deque_lock_;
	std::condition_variable log_sem_;

	static void ListenThread(GCPLogger *logger);
	std::thread log_sender_;
	volatile bool dead_;
};

GCPLogger::GCPLogger(int port, G3LogLevel default_level)
    : G3Logger(default_level), TrimFileNames(true), dead_(false)
{
	struct sockaddr_in addr;
	int yes;
	fd_ = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP);

	// Allow multiple listeners
	yes = 1;
#ifdef __linux__
	if (setsockopt(fd_, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(yes)) < 0)
		perror("Error setting SO_REUSEADDR");
#else
	if (setsockopt(fd_, SOL_SOCKET, SO_REUSEPORT, &yes, sizeof(yes)) < 0)
		perror("Error setting SO_REUSEPORT");
#endif

	addr.sin_family = AF_INET;
	addr.sin_addr.s_addr = htonl(INADDR_ANY);
	addr.sin_port = htons(port);
	if (bind(fd_, (struct sockaddr *)&addr, sizeof(addr)) < 0) {
		perror(NULL);
		dead_ = true;
		return;
	}

	if (listen(fd_, 5) < -1) {
		perror(NULL);
		dead_ = true;
		return;
	}

	log_sender_ = std::thread(ListenThread, this);
}

GCPLogger::~GCPLogger()
{
	if (dead_)
		return;

	log_deque_lock_.lock();
	dead_ = true;
	log_deque_lock_.unlock(); // Implicit memory barrier
	log_sem_.notify_all();
	log_sender_.join();

	close(fd_);
}

void
GCPLogger::Log(G3LogLevel level, const std::string &unit,
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

	std::string trimmed_filename;
	size_t lastslash = file.rfind('/');
	if (lastslash != std::string::npos && TrimFileNames)
		trimmed_filename = file.substr(lastslash+1);
	else
		trimmed_filename = file;

	int messagesize = snprintf(NULL, 0,
	    "%s (%s): %s (%s:%d in %s)",
	    log_description, unit.c_str(), message.c_str(),
	    trimmed_filename.c_str(), line, func.c_str());

	char log_message[messagesize + 1];

	sprintf(log_message,
	    "%s (%s): %s (%s:%d in %s)",
	    log_description, unit.c_str(), message.c_str(),
	    trimmed_filename.c_str(), line, func.c_str());

	log_deque_lock_.lock();
	log_deque_.push_back(log_message);

	// These are for human viewing. No one will look at > 100 messages.
	if (log_deque_.size() > 100)
		log_deque_.pop_front();
	log_sem_.notify_one();
	log_deque_lock_.unlock();
}

void GCPLogger::ListenThread(GCPLogger *logger)
{
	std::unique_lock<std::mutex> lock(logger->log_deque_lock_);
	uint64_t loglen;
	std::string curlog;

	fd_set wfds;
	int nready;
	struct timeval tv;

	do {
		lock.unlock();

		tv.tv_sec = 1;
		tv.tv_usec = 0;

		FD_ZERO(&wfds);
		FD_SET(logger->fd_, &wfds);

		nready = select(logger->fd_ + 1, NULL, &wfds, NULL, &tv);

		if (nready == 0) {
			lock.lock();
			continue;
		} else if (nready < 0) {
			fprintf(stderr, "GCPLogger select() failure: %s\n",
			    strerror(errno));
			return;
		}

		lock.lock();
		while (logger->log_deque_.empty() && !logger->dead_)
			logger->log_sem_.wait(lock);
		if (logger->dead_) {
			lock.unlock();
			return;
		}
		
		curlog = logger->log_deque_.front();
		logger->log_deque_.pop_front();
		lock.unlock();

		loglen = htonl(curlog.size());
		(void)write(logger->fd_, &loglen, sizeof(loglen));
		(void)write(logger->fd_, curlog.c_str(), curlog.size());

		lock.lock();
	} while (nready >= 0 && !logger->dead_);

	lock.unlock();
}

PYBINDINGS("gcp") {
        using namespace boost::python;

	class_<GCPLogger, bp::bases<G3Logger>, boost::shared_ptr<GCPLogger>,
	  boost::noncopyable>("GCPLogger",
	  "Logger that relays error messages to the GCP mediator over TCP",
	  init<int, G3LogLevel>((arg("port")=50030,
	       arg("default_loglevel")=G3DefaultLogLevel)))
	    .def_readwrite("trim_file_names", &GCPLogger::TrimFileNames,
	      "Show only file leaves")
	;

}

