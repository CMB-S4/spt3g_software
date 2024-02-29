#!/usr/bin/env python

import spt3g, os
from spt3g import core

p = core.G3Pipeline()
p.Add(core.G3InfiniteSource, type=core.G3FrameType.Timepoint, n=10)
p.Add(core.Dump)
p.Add(core.G3Writer, filename='testpi.g3')
p.Run()

# Check that pipelineinfo is added and is runnable
assert(len(list(core.G3File('testpi.g3'))) == 11)
for i in core.G3File('testpi.g3'):
	if i.type == core.G3FrameType.PipelineInfo:
		pi = i.values()[0]
		break

os.remove('testpi.g3')
print(repr(pi))
exec(repr(pi))
pipe.Run()

assert(len(list(core.G3File('testpi.g3'))) == 11)

# Check that PI frame has two entries on the second run through
p = core.G3Pipeline()
p.Add(core.G3Reader, filename='testpi.g3')
def check(fr):
	if fr.type == core.G3FrameType.PipelineInfo:
		assert(len(fr) == 2)
p.Add(check)
p.Run()

# Check that PI frame has two entries on the second run through even if a header
# is inserted
print('Test for twiddling')
p = core.G3Pipeline()
p.Add(core.G3Reader, filename='testpi.g3')
def twiddle(fr):
	if fr.type == core.G3FrameType.PipelineInfo:
		return [core.G3Frame(core.G3FrameType.Observation), fr]
p.Add(twiddle)
p.Add(core.G3Writer, filename='testpi2.g3')
p.Run()

# Check that it has three entries when the headered version is run from a file
print('Test for file reading after twiddling')
p = core.G3Pipeline()
p.Add(core.G3Reader, filename='testpi2.g3')
p.Add(core.Dump)
def check(fr):
	if fr.type == core.G3FrameType.PipelineInfo:
		assert(len(fr) == 3)
p.Add(check)
p.Run()

# Check that only one PipelineInfo frame is present, even if it occurs in
# multiple input files.
print('Dedup test 1')
p = core.G3Pipeline()
p.Add(core.G3Reader, filename=['testpi2.g3', 'testpi2.g3'])
p.Add(core.Dump)
i = 0
def check(fr):
	global i
	if fr.type == core.G3FrameType.PipelineInfo:
		assert(len(fr) == 3)
		i += 1
p.Add(check)
p.Run()
assert(i == 1)

# Check that deduplication does not actually deduplicate if there are
# *different* pipeline info frames in multiple input files.
print('Dedup test 2')
p = core.G3Pipeline()
p.Add(core.G3Reader, filename=['testpi.g3', 'testpi2.g3'])
p.Add(core.Dump)
i = 0
def check(fr):
	global i
	if fr.type == core.G3FrameType.PipelineInfo:
		i += 1
p.Add(check)
p.Run()
assert(i == 2)


os.remove('testpi.g3')
os.remove('testpi2.g3')
