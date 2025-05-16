#!/usr/bin/env python

from spt3g import core
import time

# File to disk
pipe = core.G3Pipeline()
pipe.Add(core.G3InfiniteSource, type=core.G3FrameType.Timepoint, n=10)
n = 0


def addinfo(fr):
    global n
    if fr.type != core.G3FrameType.Timepoint:
        return
    fr["time"] = core.G3Time(int(time.time() * core.G3Units.s))
    fr["count"] = n
    n += 1


pipe.Add(addinfo)
pipe.Add(core.Dump)
pipe.Add(core.G3Writer, filename="test.g3.xz", buffersize=1024)
pipe.Run()

# And back from disk
print("Reading")
pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename="test.g3.xz", buffersize=1024)
pipe.Add(core.Dump)
n = 0


def checkinfo(fr):
    global n
    if fr.type != core.G3FrameType.Timepoint:
        return
    if "time" not in fr:
        raise KeyError("time")
    if fr["count"] != n:
        raise ValueError("Out of order frame")
    n += 1


pipe.Add(checkinfo)
pipe.Run()

if n != 10:
    raise ValueError("Wrong number of frames (%d should be %d)" % (n, 10))
