from spt3g import core
import time

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader('/tmp/test2.g3'))
def printfr(fr):
	print(fr)
pipe.Add(printfr)

pipe.Run()
