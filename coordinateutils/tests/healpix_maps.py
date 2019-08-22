#!/usr/bin/env python

import numpy
from spt3g import coordinateutils

# Test initialization from sparse arrays

a = numpy.ones(1500)
b = numpy.arange(1500,dtype='int')
x = coordinateutils.HealpixSkyMap((b, a, 64), True, coordinateutils.MapCoordReference.Equatorial)
assert(x.nside == 64)
assert(x.npix_allocated == 1500)
assert(x[1499] == 1)
assert(x[1501] == 0)

# And dense
a = numpy.arange(49152)
x = coordinateutils.HealpixSkyMap(a)
assert(x.nside == 64)
assert(x[1533] == a[1533])

