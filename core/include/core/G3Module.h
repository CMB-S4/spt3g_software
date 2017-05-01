#ifndef _G3_MODULE_H
#define _G3_MODULE_H

#include <G3.h>
#include <G3Frame.h>
#include <G3Logging.h>

#include <deque>

class G3Module {
public:
	G3Module() {}
	virtual ~G3Module() {}

	// Process the frame in 'frame' and then put some number of frames into
	// 'out'. Cases of interest:
	// 1. The first module in the pipe (e.g. a file reader) gets NULL for
	//    frame and should fill frames into out. It should block until it
	//    has frames to send. If it does not add a frame to out, the pipe
	//    will interpret that as a signal to stop execution.
	// 2. For all other modules, frame will always be non-NULL. The usual
	//    case is to push frame onto the output queue after processing is
	//    finished. A frame can be dropped by not pushing it; it can be
	//    split into multiple frames by pushing multiple frames.
	virtual void Process(G3FramePtr frame, std::deque<G3FramePtr> &out) = 0;
};

G3_POINTERS(G3Module);

#endif
