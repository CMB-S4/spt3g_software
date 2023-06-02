#!/usr/bin/env python

from spt3g import core, maps

import numpy

m = maps.HealpixSkyMap(nside=256)

data = numpy.random.uniform(0, 1, size=m.size) > 0.5
mask = maps.G3SkyMapMask(m, data)

fr = core.G3Frame(core.G3FrameType.Map)
fr['mask'] = mask

w = core.G3Writer(filename='masktest.g3')
w(fr)
del w

f = core.G3File('masktest.g3').next()
recovereddata = numpy.asarray(f['mask'])

assert (recovereddata == data).all()

