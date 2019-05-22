#!/usr/bin/env python

from spt3g import core
from spt3g.core import G3FrameType

frametypes = [G3FrameType.Wiring, G3FrameType.Calibration, G3FrameType.Calibration, G3FrameType.Timepoint, G3FrameType.Observation, G3FrameType.Wiring, G3FrameType.Observation, G3FrameType.Calibration, G3FrameType.Scan, G3FrameType.Scan, G3FrameType.Wiring]

frames = []
for i, frametype in enumerate(frametypes):
	f = core.G3Frame(frametype)
	f['Seq'] = i
	frames.append(f)

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

pipe.Add(core.DropOrphanMetadata)

def getframes(fr):
	if fr.type != G3FrameType.EndProcessing and fr.type != G3FrameType.PipelineInfo:
		out.append(fr)
pipe.Add(getframes)

pipe.Run()

assert(len(out) == len(frames) - 3)
def dataframes(frs):
	return [f for f in frs if f.type == core.G3FrameType.Scan or f.type == core.G3FrameType.Timepoint]
assert(len(dataframes(out)) == len(dataframes(frames)))

