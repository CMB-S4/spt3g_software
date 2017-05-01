#!/usr/bin/env python

from spt3g import core
import time, sys, glob, os

# Clean up from previous run
for f in glob.glob('multitest*.g3'):
	os.remove(f)

# File to disk
pipe = core.G3Pipeline()
pipe.Add(core.G3InfiniteSource, type=core.G3FrameType.Timepoint, n=1000)
n = 0
def addinfo(fr):
	global n
	fr['time'] = core.G3Time(int(time.time()*core.G3Units.s))
	fr['count'] = n
	n += 1
pipe.Add(addinfo)
pipe.Add(core.G3MultiFileWriter, filename='multitest-%02u.g3', size_limit=20*1024)
pipe.Add(core.G3MultiFileWriter, filename=lambda frame,seq: 'multitest2-%02d.g3' % seq, size_limit=20*1024)
pipe.Add(core.G3MultiFileWriter, filename='multitest3-%02u.g3', size_limit=20*1024, divide_on=[core.G3FrameType.Timepoint])
pipe.Add(core.G3MultiFileWriter, filename='multitest4-%02u.g3', size_limit=2000000*1024, divide_on=lambda fr: fr['count'] % 200 == 0)
pipe.Run()

nfiles = len(glob.glob('multitest-*.g3'))
print('Seeing %d files, expecting 7 from format string constructor' % nfiles)
assert(nfiles == 7)

nfiles = len(glob.glob('multitest2-*.g3'))
print('Seeing %d files, expecting 7 from callback constructor' % nfiles)
assert(nfiles == 7)

nfiles = len(glob.glob('multitest3-*.g3'))
print('Seeing %d files, expecting 1000 from boundary division constructor' % nfiles)
assert(nfiles == 1000)

nfiles = len(glob.glob('multitest4-*.g3'))
print('Seeing %d files, expecting 5 from boundary callback constructor' % nfiles)
assert(nfiles == 5)

# Clean up
for f in glob.glob('multitest*.g3'):
	os.remove(f)

