#!/usr/bin/env python

import numpy
from spt3g import core, maps

# Test initialization from sparse arrays

a = numpy.arange(1500,dtype='float')
a[0] = -1
b = numpy.arange(1500,dtype='int')
x = maps.HealpixSkyMap((b, a, 64), True, False, maps.MapCoordReference.Equatorial)
assert(x.nside == 64)
assert(x.npix_allocated == 1500)
assert(x[1499] == 1499)
assert(x[1501] == 0)

# Test conversion between representations is lossless
# Along the way, test nonzero_pixels() once for each
# representation.

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

k,v = x.nonzero_pixels()
assert(len(k) == len(v) == 1500)
assert(set(v) == set(a)) # Lazy test that doesn't care about order

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

k,v = x.nonzero_pixels()
assert(len(k) == len(v) == 1500)
assert(set(v) == set(a))

x.indexedsparse = True # Dense to indexed
assert(x.nside == 64)
assert(x.npix_allocated == 1500)
for i in range(1,1500):
	assert(x[i] == i)
assert(x[1501] == 0)

k,v = x.nonzero_pixels()
assert(len(k) == len(v) == 1500)
assert(set(v) == set(a))

# Initialization from dense arrays
a = numpy.arange(49152)
x = maps.HealpixSkyMap(a)

for i in range(0, len(a)):
	assert(x[i] == i)

assert((numpy.asarray(x) == a).all())

# Conersion to ring-sparse again (trickiest, this makes sure we get all rings)
x.ringsparse = True
assert(x.npix_allocated == len(a) - 1) # First element was zero
for i in range(0, len(a)):
	assert(x[i] == i)

assert((numpy.asarray(x) == a).all())
assert(x.dense) # Should be dense again

# Conversion to/from flatsky maps
fm_stub = maps.FlatSkyMap(
        300, 300, core.G3Units.arcmin, proj=maps.MapProjection.ProjZEA
)
fm = maps.maputils.healpix_to_flatsky(x, map_stub=fm_stub)
x2 = maps.maputils.flatsky_to_healpix(fm, map_stub=x.Clone(False))

hitpix = numpy.asarray(x2) > 0
assert(numpy.allclose(numpy.asarray(x)[hitpix], numpy.asarray(x2)[hitpix]))
