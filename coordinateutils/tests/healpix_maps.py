#!/usr/bin/env python

import numpy
from spt3g import coordinateutils

# Test initialization from sparse arrays

a = numpy.arange(1500,dtype='float')
a[0] = -1
b = numpy.arange(1500,dtype='int')
x = coordinateutils.HealpixSkyMap((b, a, 64), True, coordinateutils.MapCoordReference.Equatorial)
assert(x.nside == 64)
assert(x.npix_allocated == 1500)
assert(x[1499] == 1499)
assert(x[1501] == 0)

# Test conversion between representations is lossless

x.dense = True # Indexed-sparse to dense
assert(x.npix_allocated == x.size)
assert(x[1499] == 1499)
assert(x[1501] == 0)

x.ringsparse = True # Dense to ring-sparse
assert(x.nside == 64)
assert(x.npix_allocated == 1500)
for i in range(1,1500):
	assert(x[i] == i)
assert(x[1501] == 0)

x.indexedsparse = True # Ring to indexed
assert(x.nside == 64)
assert(x.npix_allocated == 1500)
for i in range(1,1500):
	assert(x[i] == i)
assert(x[1501] == 0)

x.ringsparse = True # Indexed to ring
assert(x.nside == 64)
assert(x.npix_allocated == 1500)
assert(x[1499] == 1499)
assert(x[1501] == 0)

x.dense = True # Ring to dense
assert(x[1499] == 1499)
assert(x.npix_allocated == x.size)
assert(len({x[i] for i in range(x.size) if x[i] != 0}) == 1500)
assert(x[1501] == 0)

x.indexedsparse = True # Dense to indexed
assert(x.nside == 64)
assert(x.npix_allocated == 1500)
for i in range(1,1500):
	assert(x[i] == i)
assert(x[1501] == 0)

# Initialization from dense arrays
a = numpy.arange(49152)
x = coordinateutils.HealpixSkyMap(a)

for i in range(0, len(a)):
	assert(x[i] == i)

assert((numpy.asarray(x) == a).all())

# Conersion to ring-sparse again (trickiest, this makes sure we get all rings)
x.ringsparse = True
assert(x.npix_allocated == len(a))
for i in range(0, len(a)):
	assert(x[i] == i)

assert((numpy.asarray(x) == a).all())
assert(x.dense) # Should be dense again

