from spt3g import core
import time

pipe = core.G3Pipeline()

pipe.Add(core.G3InfiniteSource, type=core.G3FrameType.Timepoint, n=10)
def addtime(fr):
	fr['time'] = core.G3Time(int(time.time()*core.G3Units.s))
def printfr(fr):
	if fr.type != core.G3FrameType.EndProcessing:
		print(fr)
pipe.Add(addtime)
pipe.Add(printfr) # Could use core.Dump as well
pipe.Add(core.G3Writer, filename='/tmp/test2.g3')

pipe.Run()
