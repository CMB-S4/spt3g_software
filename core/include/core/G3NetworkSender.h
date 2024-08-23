#ifndef _G3_NETWORKSENDER_H
#define _G3_NETWORKSENDER_H

#include <deque>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>

#include <core/G3Module.h>

class G3NetworkSender : public G3Module {
public:
	G3NetworkSender(std::string hostname, int port, int max_queue_size,
	                int n_serializers=0);
	virtual ~G3NetworkSender();
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
	void Close(void);
	uint64_t FramesDroppedFromSerialization() const{ return n_serial_drops_; }
	uint64_t FramesDroppedFromSending() const{ return n_send_drops_; }
private:
	int fd_;
	int max_queue_size_;
	bool listening_;
	int n_serializers_;
	uint64_t n_serial_drops_;
	uint64_t n_send_drops_;

	using netbuf_type = boost::shared_ptr<std::vector<char>>;
	using buffer_future = std::shared_future<netbuf_type>;

	template<typename InputType>
	struct input_queue {
		input_queue():die(false) {}

		int max_size;
		std::mutex lock;
		std::condition_variable sem;
		std::deque<InputType> queue;

		bool die; // Protected by queue lock

		bool insert(InputType&& item, bool low_priority=false) {
			std::lock_guard<std::mutex> lg(lock);
			// Refuse to put anything into stopped queues, and if a maximum
			// queue size is set and would be exceeded, drop low priority data
			if (die || (max_size && low_priority && queue.size()>=(size_t)max_size))
				return false;
			queue.push_back(std::move(item));
			// If the queue is shared, only one worker needs to wake up, and if
			// it is private, there is only one in the first place.
			sem.notify_one();
			return true;
		}

		void stop() {
			std::lock_guard<std::mutex> lg(lock);
			die = true;
			// If the queue is shared, all workers need to wake up.
			sem.notify_all();
		}

		bool empty() const { return queue.empty(); }

		typename std::deque<InputType>::reference front() {
			return queue.front();
		}

		typename std::deque<InputType>::const_reference front() const {
			return queue.front();
		}

		void pop_front(){ queue.pop_front(); }
	};

	struct serialization_task {
		G3FramePtr input;
		std::promise<netbuf_type> output;
	};

	struct serializer_thread_data{
		serializer_thread_data(input_queue<serialization_task>& q):queue(q) {}

		std::thread thread;
		//shares a queue with other instances
		input_queue<serialization_task>& queue;
	};

	struct network_thread_data{
		std::thread thread;
		//has its own queue
		input_queue<buffer_future> queue;
		int fd;
	};

	struct output_rcord {
		G3FramePtr frame;
		buffer_future output;
	};

	input_queue<serialization_task> serialization_queue;
	std::vector<boost::shared_ptr<serializer_thread_data>> serializer_threads_;
	std::vector<boost::shared_ptr<network_thread_data>> network_threads_;
	std::deque<output_rcord> outstanding_frames;
	static void SerializeFrame(serialization_task& task);
	static void SerializeLoop(boost::shared_ptr<serializer_thread_data>);
	static void SendLoop(boost::shared_ptr<network_thread_data>);
	void StartThread(int fd);
	void StopAllThreads();
	void ReapDeadThreads();

	std::vector<std::pair<G3Frame::FrameType, buffer_future>> metadata_;

	SET_LOGGER("G3NetworkSender");
};

#endif
