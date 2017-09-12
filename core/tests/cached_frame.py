#!/usr/bin/env python

from spt3g import core
import sys

def AddInfo(frame):
    frame['Info'] = 12

# First, make sure the check doesn't trip when it isn't supposed to. 
pipe = core.G3Pipeline()
pipe.Add(core.G3InfiniteSource, type=core.G3FrameType.Timepoint, n=10)
pipe.Add(AddInfo)
pipe.Add(lambda f: 1/0)
try:
    pipe.Run() # Will throw an exception if test fails
except:
    pass

if (pipe.last_frame['Info'] != 12):
    sys.exit(1)

pipe = core.G3Pipeline()

if not pipe.last_frame is None:
    sys.exit(1)
