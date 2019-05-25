#!/usr/bin/env python

from spt3g import core
import time, sys

# File to disk
pipe = core.G3Pipeline()
pipe.Add(core.G3InfiniteSource, type=core.G3FrameType.Timepoint, n=10)
n = 0
def addinfo(fr):
	global n
	if fr.type != core.G3FrameType.Timepoint:
		return
	fr['time'] = core.G3Time(int(time.time()*core.G3Units.s))
	fr['count'] = n
	n += 1
pipe.Add(addinfo)
pipe.Add(core.Dump)
pipe.Add(core.G3Writer, filename='test.g3')
pipe.Run()

# And back from disk
print('Reading')
pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename='test.g3')
pipe.Add(core.Dump)
n = 0
def checkinfo(fr):
	global n
	if fr.type != core.G3FrameType.Timepoint:
		return
	if 'time' not in fr:
		print('No time key in frame')
		sys.exit(1)
	if fr['count'] != n:
		print('Out of order frame')
		sys.exit(1)
	n += 1
pipe.Add(checkinfo)
pipe.Run()

if n != 10:
	print('Wrong number of frames (%d should be %d)' % (n, 10))
	sys.exit(1)

sys.exit(0)

