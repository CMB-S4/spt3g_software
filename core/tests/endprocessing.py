#!/usr/bin/env python

# Test that the checks for EndProcessing frames in G3Pipeline work.

from spt3g import core
import sys

# First, make sure the check doesn't trip when it isn't supposed to.
pipe = core.G3Pipeline()
pipe.Add(core.G3InfiniteSource, type=core.G3FrameType.Timepoint, n=10)
pipe.Add(lambda f: True)
pipe.Run() # Will throw an exception if test fails

# Next, make sure it fails if EndProcessing frames are dropped
pipe = core.G3Pipeline()
pipe.Add(core.G3InfiniteSource, type=core.G3FrameType.Timepoint, n=10)
pipe.Add(lambda f: [])
try:
    pipe.Run()
except RuntimeError:
    print('Exception thrown on all frames dropped - good!')
    pass
else:
    print('No exception when EndProcessing frame dropped')
    sys.exit(1)
	
# A subtler test: make sure EndProcessing is always the *last* frame
pipe = core.G3Pipeline()
pipe.Add(core.G3InfiniteSource, type=core.G3FrameType.Timepoint, n=10)
pipe.Add(lambda f: [f, core.G3Frame(core.G3FrameType.Scan)])
try:
    pipe.Run()
except RuntimeError:
    print('Exception thrown when EndProcessing followed by scan - good!')
    pass
else:
    print('No exception when EndProcessing followed by scan')
    sys.exit(1)
	
print('EndProcessing tests pass!')
