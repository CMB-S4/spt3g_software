#include <G3EventBuilder.h>
#include <G3Pipeline.h>
#include <pybindings.h>
#include <core/SetThreadName.h>

G3EventBuilder::G3EventBuilder(int warn_size) :
  G3Module(), warn_size_(warn_size), dead_(false)
{
	process_thread_ = std::thread(ProcessThread, this);
}

G3EventBuilder::~G3EventBuilder()
{
	dead_ = true;
	process_sem_.notify_all();
	process_thread_.join();
}

void G3EventBuilder::AsyncDatum(G3TimeStamp timestamp,
    G3FrameObjectConstPtr datum)
{
	{
		std::lock_guard<std::mutex> lock(queue_lock_);
		queue_.push_back(G3EventQueueElement(timestamp, datum));
	}

	process_sem_.notify_one();
}

void G3EventBuilder::AddPolledDataModule(G3ModulePtr mod)
{
	polled_sources_.push_back(mod);
}

void G3EventBuilder::Process(G3FramePtr frame, std::deque<G3FramePtr> &out)
{
	G3PythonContext ctx("G3EventBuilder", false);

	std::unique_lock<std::mutex> lock(out_queue_lock_);

	while (out_queue_.empty() && !dead_)
		out_queue_sem_.wait(lock);

	// If dead, out_queue_ will be empty, putting nothing in the queue
	// and ending processing.

	out.swap(out_queue_);
}

void G3EventBuilder::FrameOut(G3FramePtr frame)
{
	std::lock_guard<std::mutex> lock(out_queue_lock_);

	out_queue_.push_back(frame);
	out_queue_sem_.notify_one();

	if (out_queue_.size() > 1 && out_queue_.size() % warn_size_ == 0) {
		std::string cur_mod = G3Pipeline::GetCurrentModule();
		if (cur_mod == "")
			log_warn("Outbound frame queue at %zd frames. "
			    "Possible IO stall? Rerun with profile=True to "
			    "print where.", out_queue_.size());
		else
			log_warn("Outbound frame queue at %zd frames. "
			    "Possible IO stall in module %s.",
			    out_queue_.size(), cur_mod.c_str());
	}
}

void G3EventBuilder::CollectPolledData(G3FramePtr frame)
{
	std::deque<G3FramePtr> queue, outqueue;
	queue.push_back(frame);

	for (auto i = polled_sources_.begin(); i != polled_sources_.end(); i++){
		outqueue.clear();
		for (auto j = queue.begin(); j != queue.end(); j++)
			(*i)->Process(*j, outqueue);
		queue.swap(outqueue);
	}

	if (queue.size() != 1)
		log_fatal("Need to return only 1 frame");
	if (frame != *queue.begin()) // If frame was replaced, follow changes
		*frame = **queue.begin();
}

void G3EventBuilder::ProcessThread(G3EventBuilder *builder)
{
	setThreadName("event builder");
	std::unique_lock<std::mutex> lock(builder->queue_lock_);

	while (1) {
		// If currently empty, atomically release lock and wait for data
		while (builder->queue_.empty() && !builder->dead_)
			builder->process_sem_.wait(lock);
		if (builder->dead_)
			return;

		// ProcessNewData is called with no locks, so release it
		lock.unlock();
		builder->ProcessNewData();
		lock.lock();
	}
}


PYBINDINGS("core") {
	py::class_<G3EventBuilder, py::bases<G3Module>, G3EventBuilderPtr,
	  boost::noncopyable>("G3EventBuilder", py::no_init)
	;
}

