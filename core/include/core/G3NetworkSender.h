#ifndef _G3_NETWORKSENDER_H
#define _G3_NETWORKSENDER_H

#include <deque>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>

class G3NetworkSender : public G3Module {
public:
	G3NetworkSender(std::string hostname, int port, int max_queue_size);
	virtual ~G3NetworkSender();
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
private:
	int fd_;
	int max_queue_size_;
	bool listening_;

	typedef boost::shared_ptr<std::vector<char> > netbuf_type;

	struct thread_data {
		std::thread thread;
		std::mutex queue_lock;
		std::condition_variable queue_sem;
		std::deque<netbuf_type> queue;

		int fd;
		bool die; // Protected by queue lock
	};

	std::vector<boost::shared_ptr<thread_data> > threads_;
	static void SendLoop(boost::shared_ptr<struct thread_data>);
	void StartThread(int fd);
	void ReapDeadThreads(void);

	std::vector<std::pair<G3Frame::FrameType, netbuf_type> > metadata_;

	SET_LOGGER("G3NetworkSender");
};

#endif
