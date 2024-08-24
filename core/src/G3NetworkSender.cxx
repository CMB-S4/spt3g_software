#include <pybindings.h>
#include <G3NetworkSender.h>

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <chrono>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <fcntl.h>

#include <core/SetThreadName.h>

/*
 * This class provides a unicast or broadcast distribution for frames via TCP
 * (and potentially other protocols in future). The python docstring (below
 * the class definition) provides detail on how to use it and operational
 * semantics.
 *
 * Frames are serialized to buffers either as part of the main processing
 * thread or on dedicated background threads, and then distributed to one or
 * more communication threads, each of which serves one network client. This
 * module will not output each frame until it has been completely serialized in
 * order to accomplish the following things:
 *
 * 1) Avoid paradoxes from later modification of frames by other modules.
 *    Frames, unlike frame objects, are mutable.
 * 2) Avoid worrying about races during serialization.
 *
 * Separating serialization from communication is also useful because it avoids
 * repeated serialization when serving frames to many clients.
 *
 * Using multiple serialization threads can be advantageous when frames are
 * large and high throughput is required. Otherwise, the default is to use no
 * dedicated threads for this purpose, and to simply serialize on the main
 * thread, which is sufficient in most cases.
 */


EXPORT_G3MODULE_AND("core", G3NetworkSender,
    (init<std::string, int, int, int>((arg("hostname"), arg("port"),
      arg("max_queue_size")=0, arg("n_serializers")=0))),
    "Writes frames to a network socket. If hostname is set to '*', will listen "
    "on the given port, on all interfaces, instead of connecting to the given "
    "port on a remote host. In listen mode, metadata frames (Calibration, "
    "Wiring, etc. -- everything but Scan and Timepoint) will be accumulated "
    "and the most recent of each will be sent to new clients on connect. Scan "
    "and Timepoint frames will be broadcast live to all connected clients. "
    "If max_queue_size is set to a non-zero value, Scan and Timepoint frames "
    "may be dropped if more than max_queue_size frames are queued for "
    "transmission. If n_serializers is set to a non-zero value, the task of "
    "serializing frames to be sent will be distributed across that many "
    "background threads, which is useful when high throughput of large frames "
    "is required, but is otherwise typically not necessary.",
    .def("Close", &G3NetworkSender::Close));


G3NetworkSender::G3NetworkSender(std::string hostname, int port, int max_queue,
                                 int n_serializers)
:max_queue_size_(max_queue),listening_(hostname=="*"),
 n_serializers_(n_serializers),n_serial_drops_(0),n_send_drops_(0)
{
	serialization_queue.max_size = max_queue_size_;

	if (listening_) {
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

		StartThread(fd_);
	}

	try {
		serializer_threads_.reserve(n_serializers_);
		for (size_t i=0; i<n_serializers_; i++) {
			auto t = boost::make_shared<serializer_thread_data>(serialization_queue);
			t->thread = std::thread(SerializeLoop, t);
			serializer_threads_.push_back(t);
		}
	}
	catch(...){
		//We can't rely on the destructor to do this, because if construction
		//fails, the destructor is never invoked.
		StopAllThreads();
		close(fd_);
		throw;
	}
}

G3NetworkSender::~G3NetworkSender()
{
	StopAllThreads();

	if (fd_ != -1) {
		close(fd_);
		fd_ = -1;
	}
}

void G3NetworkSender::SerializeFrame(serialization_task& task)
{
	try{
		netbuf_type buf(new std::vector<char>);
		boost::iostreams::filtering_ostream os(
			boost::iostreams::back_inserter(*buf));
		task.input->save(os);
		os.flush();
		task.output.set_value(buf);
	}
	catch (...) {
		task.output.set_exception(std::current_exception());
	}
}

void
G3NetworkSender::SerializeLoop(boost::shared_ptr<serializer_thread_data> t)
{
	setThreadName("G3NetSnd Srlize");
	std::unique_lock<std::mutex> lock(t->queue.lock);
	while(true) {
		if (t->queue.empty() && !t->queue.die)
			t->queue.sem.wait(lock);

		if (t->queue.empty()) {
			if (t->queue.die)
				break;
			//else, spurious wake-up
			continue;
		}

		serialization_task task = std::move(t->queue.front());
		t->queue.pop_front();
		lock.unlock();

		SerializeFrame(task);

		lock.lock();
	}
}

void
G3NetworkSender::SendLoop(boost::shared_ptr<network_thread_data> t)
{
	setThreadName("G3NetSnd Send");
	// Iteratively pop frames out of the queue and send them
	// on their merry way.

	std::unique_lock<std::mutex> lock(t->queue.lock);
	while (true) {
		if (t->queue.empty() && !t->queue.die)
			t->queue.sem.wait(lock);

		if (t->queue.empty()) {
			if (t->queue.die)
				break;
			//else, spurious wake-up
			continue;
		}

		std::shared_future<netbuf_type> f = t->queue.front();
		t->queue.pop_front();
		lock.unlock();

		netbuf_type buf;
		try {
			buf = f.get();
		}
		catch (...) {
			// The future might emit an exception.
			// If so, we absorb it and carry on as best we can.
			lock.lock();
			continue;
		}

		int err = write(t->fd, buf->data(), buf->size());

		lock.lock();

		// If the remote end has hung up, bail
		if (err == -1) {
			t->queue.die = true;
			break;
		}
	}
}

void
G3NetworkSender::StartThread(int fd)
{
	auto t = boost::make_shared<network_thread_data>();

	// Initialize the outbound queue with the metadata frames
	for (auto i = metadata_.begin(); i != metadata_.end(); i++)
		t->queue.queue.push_back(i->second);

	t->fd = fd;
	t->queue.max_size = max_queue_size_;
	t->thread = std::thread(SendLoop, t);

	network_threads_.push_back(t);
}

void
G3NetworkSender::StopAllThreads()
{
	serialization_queue.stop();
	for (auto& t : serializer_threads_)
		t->thread.join();
	serializer_threads_.clear();

	for (auto& t : network_threads_) {
		t->queue.stop();
		t->thread.join();
	}
	network_threads_.clear();
}

void
G3NetworkSender::ReapDeadThreads(void)
{
	// Wait for all threads marked for death to complete
	auto it=std::remove_if(network_threads_.begin(), network_threads_.end(),
		[](boost::shared_ptr<network_thread_data>& t)->bool{
			std::unique_lock<std::mutex> lock(t->queue.lock);
			if(t->queue.die){
				lock.unlock();
				t->thread.join();
				return true;
			}
			return false;
		});
	network_threads_.erase(it, network_threads_.end());
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
		StopAllThreads();

		// Wait until all frames are processed
		while (!outstanding_frames.empty()) {
			outstanding_frames.front().output.wait();
			out.push_back(outstanding_frames.front().frame);
			outstanding_frames.pop_front();
		}

		out.push_back(frame);
		return;
	} else {
		// Send other frames normally

		// Arrange getting the frame serialized
		bool non_meta = (frame->type == G3Frame::Scan ||
		                 frame->type == G3Frame::Timepoint);
		serialization_task s_task;
		s_task.input = frame;
		buffer_future s_future = s_task.output.get_future().share();
		bool will_serialize=true;
		if (n_serializers_) {
			will_serialize = serialization_queue.insert(std::move(s_task),
			                                            non_meta);
		}
		else
			SerializeFrame(s_task); // do directly on the main thread

		// Clean up after any threads that have stopped on their own
		ReapDeadThreads();

		// Each network thread gets its own copy of the future
		if (will_serialize) {
			for (auto& t : network_threads_) {
				if (not t->queue.insert(buffer_future(s_future), non_meta))
					n_send_drops_++;
			}
			outstanding_frames.push_back(output_rcord{frame,
			                                          buffer_future(s_future)});
		}
		else
			n_serial_drops_++;

		// If this is a metadata frame (not scan or timepoint),
		// stick it in the queue.
		if (!non_meta) {
			auto i = metadata_.begin();
			while (i != metadata_.end()) {
				if (i->first == frame->type) {
					i->second = s_future;
					break;
				}
				i++;
			}

			if (i == metadata_.end())
				metadata_.push_back(std::make_pair(frame->type, s_future));
		}
	}

	// Output any frames which have been processed
	while (!outstanding_frames.empty() &&
	       outstanding_frames.front().output.wait_for(std::chrono::seconds(0))
	       == std::future_status::ready) {
		out.push_back(outstanding_frames.front().frame);
		outstanding_frames.pop_front();
	}
}

void
G3NetworkSender::Close()
{
	// Ask all threads to please die now.
	StopAllThreads();

	if (listening_) {
		close(fd_);
		fd_ = -1;
	}
}
