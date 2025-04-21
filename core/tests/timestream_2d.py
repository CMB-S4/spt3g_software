#!/usr/bin/env python

import numpy
from spt3g import core

# Check that timestream maps being cast with numpy.asarray have the right
# content and are indistinguishable from numpy.asarray(list(tsm.values())

tsm = core.G3TimestreamMap()
start = core.G3Time.Now()
stop = start + 5*core.G3Units.s
keys = ['D', 'A', 'B', 'C']
for ts in keys:
	i = core.G3Timestream(numpy.random.normal(size=600))
	i.start = start
	i.stop = stop
	i.units = core.G3TimestreamUnits.Tcmb
	tsm[ts] = i

buffer1d = numpy.asarray(list(tsm.values()))

buffer2d = numpy.asarray(tsm)

numpy.testing.assert_array_equal(keys, list(tsm.keys()))
assert(buffer2d.shape == (4,600))
numpy.testing.assert_array_equal(buffer1d, buffer2d)

# Try building the TSM directly
tsm2 = core.G3TimestreamMap(tsm.keys(), buffer1d)
numpy.testing.assert_array_equal(numpy.asarray(tsm2), buffer1d)
numpy.testing.assert_array_equal(numpy.asarray(tsm2), buffer2d)

# Try building semi-directly, using a list of numpy arrays
tsm2 = core.G3TimestreamMap(tsm.keys(), [x for x in buffer1d])
numpy.testing.assert_array_equal(numpy.asarray(tsm2), buffer1d)
numpy.testing.assert_array_equal(numpy.asarray(tsm2), buffer2d)

# Round-trip to and from numpy
tsm2 = core.G3TimestreamMap(tsm.keys(), buffer2d)
numpy.testing.assert_array_equal(numpy.asarray(tsm2), buffer1d)
numpy.testing.assert_array_equal(numpy.asarray(tsm2), buffer2d)

# Now try writing to it
numpy.asarray(tsm2)[1,5] = 16
assert(numpy.asarray(tsm2)[1,5] == 16)

# And that writes propagate back to the underlying G3Timestream
buffer2d[1] = numpy.random.normal(size=600)
numpy.testing.assert_array_equal(numpy.asarray(tsm['B']), buffer2d[keys.index('B')])

# Test initialization from numpy and data copying/sharing
buffer1d[2,6] = -3.0

# By default (copy_data = True), changes should not propagate from the array used to initialize
# to the timestreams
tsm2 = core.G3TimestreamMap(tsm.keys(), buffer1d)
buffer1d[2,6] = 14.0
assert(numpy.asarray(tsm2)[2,6] != 14.0)

# But we should be able to set things up so that they do
tsm2 = core.G3TimestreamMap(tsm.keys(), buffer1d, copy_data = False)
buffer1d[2,6] = 16.0
assert(numpy.asarray(tsm2)[2,6] == 16.0)

# Test setting secondary properties at the time of direct construction
tsm2 = core.G3TimestreamMap(tsm.keys(), buffer1d, start=start, stop=stop,
                            units=core.G3TimestreamUnits.Counts,
                            compression_level=3, bit_depth=24)
for ts in tsm2.values():
	assert ts.start == start
	assert ts.stop == stop
	assert ts.units == core.G3TimestreamUnits.Counts
	assert ts.compression_level == 3
	assert ts.bit_depth == 24
