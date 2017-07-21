from spt3g import core
import time, sys

# Read in one or more files specified at the command line and print their
# contents to the console.

pipe = core.G3Pipeline()

pipe.Add(core.G3Reader, filename=sys.argv[1:])
def printfr(fr):
	print(fr)
pipe.Add(printfr) # core.Dump does the same thing as this

pipe.Run()
