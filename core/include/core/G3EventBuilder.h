#ifndef _G3_EVENTBUILDER_H
#define _G3_EVENTBUILDER_H

#include <G3Frame.h>
#include <G3Module.h>
#include <G3TimeStamp.h>

#include <deque>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>

class G3EventBuilder : public G3Module {
public:
	G3EventBuilder(int warn_size = 1000); // Yell if the output frame queue
	                                      // grows larger than warn_size
	virtual ~G3EventBuilder();

	void AsyncDatum(G3TimeStamp, G3FrameObjectConstPtr);
	void AddPolledDataModule(G3ModulePtr);

	void Process(G3FramePtr, std::deque<G3FramePtr> &);

protected:
	void CollectPolledData(G3FramePtr frame);
	virtual void ProcessNewData() = 0;
	void FrameOut(G3FramePtr frame);

	typedef std::pair<G3TimeStamp, G3FrameObjectConstPtr>
	    G3EventQueueElement;
	typedef std::deque<G3EventQueueElement> G3EventQueue;

	std::mutex queue_lock_;
	G3EventQueue queue_;

private:
	static void ProcessThread(G3EventBuilder *);

	int warn_size_;

	std::thread process_thread_;
	std::condition_variable process_sem_;
	std::vector<G3ModulePtr> polled_sources_;
	bool dead_;

	std::condition_variable out_queue_sem_;
	std::mutex out_queue_lock_;
	std::deque<G3FramePtr> out_queue_;

	SET_LOGGER("G3EventBuilder");
};

G3_POINTERS(G3EventBuilder);

#endif
