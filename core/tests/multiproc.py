#!/usr/bin/env python

from spt3g import core
import time, sys

# File to disk
n = 0
def addinfo(fr):
	global n
	if fr.type != core.G3FrameType.Timepoint:
		return
	fr['time'] = core.G3Time(int(time.time()*core.G3Units.s))
	fr['count'] = n
	n += 1
m = 0
def checkinfo(fr):
	global m
	if fr.type != core.G3FrameType.Timepoint:
		return
	if 'time' not in fr:
		print('No time key in frame')
		sys.exit(1)
	if fr['count'] != m:
		print('Out of order frame')
		sys.exit(1)
	m += 1

if __name__ == '__main__':
	pipe = core.G3Pipeline()
	pipe.Add(core.G3InfiniteSource, type=core.G3FrameType.Timepoint, n=10)
	pipe.Add(addinfo, subprocess=True)
	pipe.Add(core.Dump)
	pipe.Add(checkinfo)
	pipe.Run()

	if n is not 0:
		print('Ran in same process!')
		sys.exit(1)
	n = 10 # Different process, so don't get this back
	if m != n:
		print('Wrong number of frames (%d should be %d)' % (m, n))
		sys.exit(1)

	sys.exit(0)

