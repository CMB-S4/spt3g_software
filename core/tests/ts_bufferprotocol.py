#!/usr/bin/env python

import numpy, sys
from spt3g import core

ts = core.G3Timestream(numpy.ones(1200))

tsarr = numpy.asarray(ts)

if tsarr[5] != ts[5]:
	print('Values corrupted')
	sys.exit(1)

if len(tsarr) != len(ts):
	print('Lengths do not match')
	sys.exit(1)

tsarr[6] = 15
if ts[6] != 15:
	print('Buffer not shared')
	sys.exit(1)

ts[7] = 45
if tsarr[7] != 45:
	print('Buffer not shared')
	sys.exit(1)

print(tsarr[5:9])

sys.exit(0)

