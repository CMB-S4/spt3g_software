#!/usr/bin/env python

import pickle, sys
from spt3g import core

f = core.G3Frame(core.G3FrameType.Scan)
f['Five'] = 5

b = pickle.dumps(f)
f2 = pickle.loads(b)

assert(f2['Five'] == 5)
assert(f2.type == f.type)

sys.exit(0)

