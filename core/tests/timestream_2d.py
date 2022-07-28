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

# Try building the TSM directly
tsm2 = core.G3TimestreamMap(tsm.keys(), buffer1d)
assert((numpy.asarray(tsm2) == buffer1d).all())
assert((numpy.asarray(tsm2) == buffer2d).all())

# Try building semi-directly, using a list of numpy arrays
tsm2 = core.G3TimestreamMap(tsm.keys(), [x for x in buffer1d])
assert((numpy.asarray(tsm2) == buffer1d).all())
assert((numpy.asarray(tsm2) == buffer2d).all())

# Round-trip to and from numpy
tsm2 = core.G3TimestreamMap(tsm.keys(), buffer2d)
assert((numpy.asarray(tsm2) == buffer1d).all())
assert((numpy.asarray(tsm2) == buffer2d).all())

# Now try writing to it
buffer2d[1] = numpy.random.normal(size=600)
assert((numpy.asarray(tsm['B']) == buffer2d[1]).all())


