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
		raise RuntimeError('No time key in frame')
	if fr['count'] != m:
		raise RuntimeError('Out of order frame')
	m += 1

if __name__ == '__main__':
	if sys.version_info[:2] > (3, 7):
		print('Subprocess option is disabled for python versions > 3.7')
		raise SystemExit

	pipe = core.G3Pipeline()
	pipe.Add(core.G3InfiniteSource, type=core.G3FrameType.Timepoint, n=10)
	pipe.Add(addinfo, subprocess=True)
	pipe.Add(core.Dump)
	pipe.Add(checkinfo)
	pipe.Run()

	if n != 0:
		raise RuntimeError('Ran in same process!')
	n = 10 # Different process, so don't get this back
	if m != n:
		raise RuntimeError('Wrong number of frames (%d should be %d)' % (m, n))
