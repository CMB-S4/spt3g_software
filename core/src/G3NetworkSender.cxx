#include <pybindings.h>
#include <G3Module.h>
#include <G3NetworkSender.h>

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <deque>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <fcntl.h>

/*
 * This class provides a unicast or broadcast distribution for frames via TCP
 * (and potentially other protocols in future). The python docstring (below
 * the class definition) provides detail on how to use it and operational
 * semantics.
 *
 * Frames are serialized to buffers as part of the main processing thread
 * and are then distributed to one or more communication threads, each of which
 * serves one network client. Data are communicated to the communication threads
 * as buffers rather than frames to accomplish the following things:
 *
 * 1) Avoids paradoxes from later modification of frames by other modules.
 *    Frames, unlike frame objects, are mutable.
 * 2) Avoids worrying about races during serialization.
 * 3) Centralizes the serialization CPU load in the event that many clients
 *    are listening.
 * 4) Works around any potential need to acquire the GIL when destroying shared
 *    pointers that may have come from Python
 */


EXPORT_G3MODULE("core", G3NetworkSender,
    (init<std::string, int, int>((arg("hostname"), arg("port"),
      arg("max_queue_size")=0))),
    "Writes frames to a network socket. If hostname is set to '*', will listen "
    "on the given port, on all interfaces, instead of connecting to the given "
    "port on a remote host. In listen mode, metadata frames (Calibration, "
    "Wiring, etc. -- everything but Scan and Timepoint) will be accumulated "
    "and the most recent of each will be sent to new clients on connect. Scan "
    "and Timepoint frames will be broadcasted live to all connected clients. "
    "If max_queue_size is set to a non-zero value, Scan and Timepoint frames "
    "may be dropped if more than max_queue_size frames are queued for "
    "transmission.");

G3NetworkSender::G3NetworkSender(std::string hostname, int port, int max_queue)
{
	max_queue_size_ = max_queue;

	if (strcmp(hostname.c_str(), "*") == 0) {
		// Listen for incoming connections

		struct sockaddr_in6 sin;
		int no = 0, yes = 1;

		// Use both IPv6 and IPv4
		bzero(&sin, sizeof(sin));
		sin.sin6_family = AF_INET6;
		sin.sin6_port = htons(port);
	#ifdef SIN6_LEN
		sin.sin6_len = sizeof(sin);
	#endif

		fd_ = socket(PF_INET6, SOCK_STREAM, 0);
		if (fd_ <= 0)
			log_fatal("Could not listen on port %d (%s)",
			    port, strerror(errno));
		setsockopt(fd_, IPPROTO_IPV6, IPV6_V6ONLY, &no,
		    sizeof(no));
		setsockopt(fd_, SOL_SOCKET, SO_REUSEADDR, &yes,
		    sizeof(yes));

		// Listening socket needs to be nonblocking so that we can
		// check if connections are waiting.
		fcntl(fd_, F_SETFL, fcntl(fd_, F_GETFL, 0) | O_NONBLOCK);

		if (bind(fd_, (struct sockaddr *)&sin, sizeof(sin)) < 0)
			log_fatal("Could not bind on port %d (%s)",
			    port, strerror(errno));
		if (listen(fd_, 10) < 0)
			log_fatal("Could not listen on port %d (%s)",
			    port, strerror(errno));

		listening_ = true;
	} else {
		// Connect to a listening host elsewhere

		struct addrinfo hints, *info, *r;
		char portstr[16];
		int err;

		bzero(&hints, sizeof(hints));
		hints.ai_family = AF_UNSPEC;
		hints.ai_socktype = SOCK_STREAM;
		snprintf(portstr, sizeof(portstr), "%d", port);

		err = getaddrinfo(hostname.c_str(), portstr, &hints, &info);
		if (err != 0)
			log_fatal("Could not find host %s (%s)",
			    hostname.c_str(), gai_strerror(err));

		// Loop through possible addresses until we find one
		// that works.
		fd_ = -1;
		for (r = info; r != NULL; r = r->ai_next) {
			fd_ = socket(r->ai_family, r->ai_socktype,
			    r->ai_protocol);
			if (fd_ == -1)
				continue;

			if (connect(fd_, r->ai_addr, r->ai_addrlen) == -1) {
				close(fd_);
				fd_ = -1;
				continue;
			}

			break;
		}

		if (fd_ == -1)
			log_fatal("Could not connect to %s:%d (%s)",
			    hostname.c_str(), port, strerror(errno));

		if (info != NULL)
			freeaddrinfo(info);

		listening_ = false;
		StartThread(fd_);
	}
}

G3NetworkSender::~G3NetworkSender()
{

	for (auto i = threads_.begin(); i != threads_.end(); i++) {
		(*i)->queue_lock.lock();
		(*i)->die = true;
		(*i)->queue_sem.notify_one();
		(*i)->queue_lock.unlock();

		(*i)->thread.join();
	}
}

void
G3NetworkSender::SendLoop(boost::shared_ptr<struct thread_data> t)
{
	// Iteratively pop frames out of the queue and send them
	// on their merry way.

	std::unique_lock<std::mutex> lock(t->queue_lock);
	int err;

	while (true) {
		if (t->queue.empty() && !t->die)
			t->queue_sem.wait(lock);

		if (t->queue.empty() && t->die)
			break;

		netbuf_type buf = t->queue.front();
		t->queue.pop_front();
		lock.unlock();

		err = write(t->fd, buf->data(), buf->size());

		lock.lock();

		// If the remote end has hung up, bail
		if (err == -1) {
			t->die = true;
			break;
		}
	}
}

void
G3NetworkSender::StartThread(int fd)
{
	boost::shared_ptr<thread_data> t(new thread_data);

	// Initialize the outbound queue with the metadata frames
	for (auto i = metadata_.begin(); i != metadata_.end(); i++)
		t->queue.push_back(i->second);

	t->fd = fd;
	t->die = false;
	t->thread = std::thread(SendLoop, t);

	threads_.push_back(t);
}

void
G3NetworkSender::ReapDeadThreads(void)
{
	// Wait for all threads marked for death to complete
restart:
	for (auto i = threads_.begin(); i != threads_.end(); i++) {
		(*i)->queue_lock.lock();
		if ((*i)->die) {
			(*i)->queue_lock.unlock();
			(*i)->thread.join();
			threads_.erase(i);
			goto restart;
		}

		(*i)->queue_lock.unlock();
	}
}


void
G3NetworkSender::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
{
	if (listening_) {
		// Check if we have any pending connections
		int fd;

		while ((fd = accept(fd_, NULL, NULL)) != -1) {
			// Actual transfer socket needs to be set back to
			// blocking
			fcntl(fd, F_SETFL, fcntl(fd, F_GETFL, 0) & ~O_NONBLOCK);
			StartThread(fd);
		}
	}

	if (frame->type == G3Frame::EndProcessing) {
		// Reap all threads on EndProcessing.
		for (auto i = threads_.begin(); i != threads_.end(); i++) {
			(*i)->queue_lock.lock();
			(*i)->die = true;
			(*i)->queue_sem.notify_one();
			(*i)->queue_lock.unlock();
		}
		ReapDeadThreads();
	} else {
		// Send other frames normally
		netbuf_type buf(new std::vector<char>);
		{
			boost::iostreams::filtering_ostream os(
			    boost::iostreams::back_inserter(*buf));
			frame->save(os);
			os.flush();
		}

		// Clean up after any threads that have stopped on their own
		ReapDeadThreads();

		for (auto i = threads_.begin(); i != threads_.end(); i++) {
			(*i)->queue_lock.lock();

			// Check if we need to drop some packets.
			if (max_queue_size_ > 0 &&
			   (*i)->queue.size() > max_queue_size_ &&
			    (frame->type == G3Frame::Scan ||
			     frame->type == G3Frame::Timepoint)) {
				(*i)->queue_lock.unlock();
				continue;
			}

			// Otherwise, send away
			(*i)->queue.push_back(buf);
			(*i)->queue_lock.unlock();

			(*i)->queue_sem.notify_one();
		}

		// If this is a metadata frame (not scan or timepoint),
		// stick it in the queue.
		if (frame->type != G3Frame::Scan &&
		    frame->type != G3Frame::Timepoint) {
			auto i = metadata_.begin();
			while (i != metadata_.end()) {
				if (i->first == frame->type) {
					i->second = buf;
					break;
				}
				i++;
			}

			if (i == metadata_.end())
				metadata_.push_back(std::pair<
				    G3Frame::FrameType, netbuf_type>(
				    frame->type, buf));
		}
	}
	out.push_back(frame);
}
