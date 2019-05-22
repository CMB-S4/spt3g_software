#!/usr/bin/env python

from spt3g import core
import sys

# File to disk
pipe = core.G3Pipeline()
pipe.Add(core.G3InfiniteSource, type=core.G3FrameType.Timepoint, n=10)
n = 0
def addinfo(fr):
	global n
	if fr.type != core.G3FrameType.Timepoint:
		return
	fr['count'] = n
	n += 1
pipe.Add(addinfo)
pipe.Add(lambda frame: 'count' in frame and frame['count'] > 5)
pipe.Add(core.Dump)
ids = []
def getids(fr):
	if fr.type != core.G3FrameType.Timepoint:
		return
	ids.append(fr['count'])
pipe.Add(getids)

pipe.Run()

print(ids)
if ids != [6,7,8,9]:
	sys.exit(1)

sys.exit(0)

