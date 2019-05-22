#!/usr/bin/env python

from spt3g import core

frames = []
frames.append(core.G3Frame(core.G3FrameType.Calibration))
frames.append(core.G3Frame(core.G3FrameType.Calibration))
frames.append(core.G3Frame(core.G3FrameType.Scan))
frames.append(core.G3Frame(core.G3FrameType.Calibration))
frames.append(core.G3Frame(core.G3FrameType.Calibration))
frames.append(core.G3Frame(core.G3FrameType.Scan))
frames.append(core.G3Frame(core.G3FrameType.Scan))

# Make first two cal frames different, the rest identical
frames[0]['Seq'] = 'B'
frames[1]['Seq'] = 'A'
frames[3]['Seq'] = 'A'
frames[4]['Seq'] = 'A'

# Make all Scan frames identical
frames[2]['Seq'] = 'B'
frames[5]['Seq'] = 'B'
frames[6]['Seq'] = 'B'

# Result should be that we get the first two calibration frames, since they are
# different, all the scan frames, since it should not even try to deduplicate
# data, and nothing else

out = []
pipe = core.G3Pipeline()

sent = False
def framesource(fr):
	global sent
	if not sent:
		sent = True
		return frames
	else:
		return []
pipe.Add(framesource)

pipe.Add(core.DeduplicateMetadata)

def getframes(fr):
	if fr.type != core.G3FrameType.EndProcessing and fr.type != core.G3FrameType.PipelineInfo:
		out.append(fr)
pipe.Add(getframes)

pipe.Run()

assert(len(out) == 5)
assert(out[0].type == core.G3FrameType.Calibration)
assert(out[0]['Seq'] == 'B')
assert(out[1].type == core.G3FrameType.Calibration)
assert(out[1]['Seq'] == 'A')


assert(out[2].type == core.G3FrameType.Scan)
assert(out[3].type == core.G3FrameType.Scan)
assert(out[4].type == core.G3FrameType.Scan)

