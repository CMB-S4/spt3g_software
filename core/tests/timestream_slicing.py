#!/usr/bin/env python

import numpy
from spt3g import core

ts = core.G3Timestream(numpy.zeros(50))
ts.units = core.G3TimestreamUnits.Power
ts.start = core.G3Time(0)
ts.stop = core.G3Time(10*core.G3Units.s)

assert(ts[5] == 0) # Check scalar getitem

ts[12] = 13.4 # Check scalar setitem
assert(ts[12] == 13.4)

ts[:] = numpy.arange(50) # Test vector setitem
assert(ts[12]) == 12.0

# Test units preserved by slicing
assert(ts[:5].units == ts.units)

# Test length with slicing
assert(len(ts[::2]) == len(ts)/2)
assert(len(ts[5:]) == len(ts)-5)

# Test consistency with numpy slicing with steps that do and do not divide evenly into array length
assert((numpy.asarray(ts[::2]) == numpy.asarray(ts)[::2]).all())
assert((numpy.asarray(ts[::3]) == numpy.asarray(ts)[::3]).all())
assert((numpy.asarray(ts[5::2]) == numpy.asarray(ts)[5::2]).all())
assert((numpy.asarray(ts[5::3]) == numpy.asarray(ts)[5::3]).all())
assert((numpy.asarray(ts[5:-4:2]) == numpy.asarray(ts)[5:-4:2]).all())
assert((numpy.asarray(ts[:-4:2]) == numpy.asarray(ts)[:-4:2]).all())
assert((numpy.asarray(ts[:-4]) == numpy.asarray(ts)[:-4]).all())
assert((numpy.asarray(ts[5:]) == numpy.asarray(ts)[5:]).all())
assert((numpy.asarray(ts[5:-4]) == numpy.asarray(ts)[5:-4]).all())

# Test sample rate consistency for the first and last N samples
assert(numpy.isclose(ts[:5].sample_rate, ts.sample_rate, atol=0))
assert(numpy.isclose(ts[-5:].sample_rate, ts.sample_rate, atol=0))

# Test start and stop time accuracy
assert(numpy.allclose([t.time for t in ts[:5].times()], [t.time for t in ts.times()][:5]))
assert(numpy.allclose([t.time for t in ts[-10:].times()], [t.time for t in ts.times()][-10:]))
assert(numpy.allclose([t.time for t in ts[5:-10].times()], [t.time for t in ts.times()][5:-10]))

# Test that we get the right sample rates when sparsely sampling
assert(numpy.isclose(ts[::5].sample_rate, ts.sample_rate/5, atol=0))
assert(numpy.isclose(ts[::2].sample_rate, ts.sample_rate/2, atol=0))

# Mix all the slicing tests together in a big bucket
assert(numpy.isclose(ts[12::2].sample_rate, ts.sample_rate/2, atol=0))
assert(numpy.isclose(ts[12:-2:2].sample_rate, ts.sample_rate/2, atol=0))
assert(numpy.isclose(ts[:-2:2].sample_rate, ts.sample_rate/2, atol=0))

