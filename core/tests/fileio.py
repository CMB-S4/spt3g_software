#!/usr/bin/env python

from spt3g import core
import time

nframes = 1000
verbose = False

# File to disk
core.log_notice("Writing %d frames to disk" % nframes)
pipe = core.G3Pipeline()
pipe.Add(core.G3InfiniteSource, type=core.G3FrameType.Timepoint, n=nframes)
# Drop a wiring frame in the middle
count = 0
def addwiring(fr):
    global count
    if fr.type == core.G3FrameType.Timepoint:
        count += 1
        if count == 6:
            fr2 = core.G3Frame(core.G3FrameType.Wiring)
            return [fr2, fr]
pipe.Add(addwiring)
n = 0
def addinfo(fr):
    global n
    if fr.type != core.G3FrameType.Timepoint:
        return
    fr['time'] = core.G3Time(int(time.time()*core.G3Units.s))
    fr['count'] = n
    n += 1
pipe.Add(addinfo)
if verbose:
    pipe.Add(core.Dump)
pipe.Add(core.G3Writer, filename='test.g3')
pipe.Run(profile=True)
assert n == nframes, 'Wrong number of frames written (%d should be %d)' % (n, nframes)

# And back from disk
core.log_notice("Reading %d frames from disk" % nframes)
pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename='test.g3', track_filename=True)
if verbose:
    pipe.Add(core.Dump)
n = 0
def checkinfo(fr):
    global n
    if fr.type != core.G3FrameType.Timepoint:
        return
    assert 'time' in fr, 'No time key in frame'
    assert fr['count'] == n, 'Out of order frame'
    assert fr._filename == 'test.g3', 'Wrong filename'
    n += 1
pipe.Add(checkinfo)
pipe.Run(profile=True)

assert n == nframes, 'Wrong number of frames read (%d should be %d)' % (n, nframes)

# Skip empty files
core.log_notice("Reading %d frames from disk with empty files" % nframes)
with core.G3Writer("empty.g3") as wr:
    pass
n = 0
pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename=["empty.g3", "test.g3", "empty.g3"], track_filename=True)
pipe.Add(checkinfo)
pipe.Run(profile=True)

assert n == nframes, 'Wrong number of frames read (%d should be %d)' % (n, nframes)

# Indexing
class CachingReader:
    def __init__(self, filename='test.g3'):
        self.reader = core.G3Reader(filename='test.g3')
        self.w_pos = None

    def __call__(self, frame):
        assert frame is None
        pos = self.reader.tell()
        fr = self.reader(frame)
        if not len(fr):
            return fr
        if fr[0].type == core.G3FrameType.Wiring:
            self.w_pos = pos
        return fr

cacher = CachingReader()

core.log_notice("Reading %d frames from disk with frame indexing" % nframes)
pipe = core.G3Pipeline()
pipe.Add(cacher)
if verbose:
    pipe.Add(core.Dump)
pipe.Run(profile=True)

assert cacher.w_pos is not None, 'Missing wiring frame'

# Using cached index
class CachedReader:
    def __init__(self, filename='test.g3', start=None):
        self.reader = core.G3Reader(filename=filename)
        self.pos = start

    def __call__(self, frame):
        assert frame is None
        if self.pos is not None:
            self.reader.seek(self.pos)
        if self.pos is not None:
            assert self.reader.tell() == self.pos
        fr = self.reader(frame)
        if not len(fr):
            return fr
        if self.pos is not None:
            assert fr[0].type == core.G3FrameType.Wiring
            self.pos = None
        return fr

cached = CachedReader(start=cacher.w_pos)

core.log_notice("Reading %d frames from disk with seek" % nframes)
pipe = core.G3Pipeline()
pipe.Add(cached)
if verbose:
    pipe.Add(core.Dump)
pipe.Run(profile=True)

assert cached.pos is None, "Missing wiring frame"

# check EOF handling
while cached.reader(None):
    pass

# can seek to EOF if at EOF
pos = cached.reader.tell()
cached.reader.seek(pos)

# can't seek backward once at EOF
try:
    cached.reader.seek(0)
except:
    pass
else:
    assert False, "Expected RuntimeError seeking back from EOF"
