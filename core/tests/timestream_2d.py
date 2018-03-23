#!/usr/bin/env python

import numpy
from spt3g import core

# Check that timestream maps being cast with numpy.asarray have the right
# content and are indistinguishable from numpy.asarray(list(tsm.values())

tsm = core.G3TimestreamMap()
start = core.G3Time.Now()
stop = start + 5*core.G3Units.s
for ts in ['A', 'B', 'C', 'D']:
	i = core.G3Timestream(numpy.random.normal(size=600))
	i.start = start
	i.stop = stop
	i.units = core.G3TimestreamUnits.Tcmb
	tsm[ts] = i

buffer1d = numpy.asarray(list(tsm.values()))

buffer2d = numpy.asarray(tsm)

assert(buffer1d.shape == buffer2d.shape)
assert(buffer2d.shape == (4,600))
assert((buffer1d == buffer2d).all())

