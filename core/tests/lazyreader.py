#!/usr/bin/env python
from spt3g import core
import sys, glob

class NumberAdder(object):
    def __init__(self):
        self.n = 0
    def __call__(self, frame):
        if frame.type == core.G3FrameType.EndProcessing:
            return
        frame['TheNumber'] = self.n
        self.n += 1

#set up our file to read
pipe = core.G3Pipeline()
pipe.Add(core.G3InfiniteSource, type=core.G3FrameType.Timepoint, n=10)
pipe.Add(core.G3Writer, filename = 'testlazyreader.g3')
pipe.Run()

na = NumberAdder()
pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename = 'testlazyreader.g3', n_frames_to_read = 3)
pipe.Add(na)
pipe.Run() 

if na.n != 3:
    print("N Reader failed", na.n)
    sys.exit(1)


