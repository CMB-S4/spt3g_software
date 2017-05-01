import sys
from spt3g import core

pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename = sys.argv[1])
pipe.Add(core.G3Writer, filename = sys.argv[2])
pipe.Run()
