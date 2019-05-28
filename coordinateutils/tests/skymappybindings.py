#!/usr/bin/env python

import numpy
from spt3g.coordinateutils import G3SkyMap, MapCoordReference

# 1D arrays
m = G3SkyMap(MapCoordReference.Equatorial, xpixels=500)
assert(len(m.shape) == 1)
assert(m.shape[0] == 500)
m[15] = 65.4
assert(numpy.asarray(m)[15] == 65.4)
assert(numpy.asarray(m).shape == m.shape)

# 2D arrays
m = G3SkyMap(MapCoordReference.Equatorial, xpixels=500, ypixels=20)
assert(len(m.shape) == 2)
assert(m.shape[0] == 20)
assert(m.shape[1] == 500)
assert(m.shape[0] == numpy.asarray(m).shape[0])
assert(m.shape[1] == numpy.asarray(m).shape[1])
m[15] = 65.4
assert(numpy.asarray(m)[0,15] == 65.4)
m[3,15] = 65.4
assert(numpy.asarray(m)[3,15] == 65.4)
assert(numpy.asarray(m).shape == m.shape)

# Reverse direction 2D arrays
w = numpy.random.uniform(size=(300,50))
m = G3SkyMap(w, MapCoordReference.Equatorial)
assert(len(m.shape) == 2)
assert(m.shape[0] == w.shape[0])
assert(m.shape[1] == w.shape[1])
assert(w[17,32] == m[17,32])
assert((w == m).all())

