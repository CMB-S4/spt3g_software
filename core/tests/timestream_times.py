#!/usr/bin/env python

from __future__ import print_function

import numpy, sys
from spt3g import core

then = core.G3Time.Now()
now = core.G3Time(then.time + 3*core.G3Units.second)

data = numpy.zeros(200)

ts = core.G3Timestream(data)
ts.start = then
ts.stop = now

assert(numpy.abs(ts.sample_rate/core.G3Units.Hz - 66.33333) < 1e-5)

times = ts.times()

assert(ts.times()[0] == ts.start)
print(ts.times()[-1].time, ts.stop.time)
assert(numpy.abs(ts.times()[-1].time - ts.stop.time) < 2) # Max 1 tick roundoff error

tsm = core.G3TimestreamMap()
tsm['Test'] = ts

assert(numpy.abs(tsm.sample_rate/core.G3Units.Hz - 66.33333) < 1e-5)

times = tsm.times()

assert(tsm.times()[0] == tsm.start)
assert(numpy.abs(tsm.times()[-1].time - tsm.stop.time) < 2) # Max 1 tick roundoff error

