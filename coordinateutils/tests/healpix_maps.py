#!/usr/bin/env python

import numpy
from spt3g import core, coordinateutils

# Test initialization from sparse arrays

a = numpy.arange(1500,dtype='float')
a[0] = -1
b = numpy.arange(1500,dtype='int')
x = coordinateutils.HealpixSkyMap((b, a, 64), True, False, coordinateutils.MapCoordReference.Equatorial)
assert(x.nside == 64)
assert(x.npix_allocated == 1500)
assert(x[1499] == 1499)
assert(x[1501] == 0)

# Test conversion between representations is lossless
# Along the way, test nonzero_pixels() and rebin() once for each
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

import healpy as hp
p0 = hp.nest2ring(64, (hp.ring2nest(32, 0) * 4 + numpy.arange(4)).astype(int))
v0 = sum([x[int(i)] for i in p0])

x2 = x.rebin(2, norm=False)
assert(x2[0] == v0)
assert(numpy.sum(x2) == numpy.sum(v))

x.shift_ra = True
x.ringsparse = True # Indexed to ring
assert(x.nside == 64)
assert(x.npix_allocated == 1512) # shifted ringsparse is less efficient
assert(x[1499] == 1499)
assert(x[1501] == 0)

k,v = x.nonzero_pixels()
assert(len(k) == len(v) == 1500)
assert(set(v) == set(a))

x2 = x.rebin(2, norm=False)
assert(x2[0] == v0)
assert(numpy.sum(x2) == numpy.sum(v))

x.dense = True # Ring to dense
x.shift_ra = False
assert(x[1499] == 1499)
assert(x.npix_allocated == x.size)
assert(len({x[i] for i in range(x.size) if x[i] != 0}) == 1500)
assert(x[1501] == 0)

k,v = x.nonzero_pixels()
assert(len(k) == len(v) == 1500)
assert(set(v) == set(a))

x2 = x.rebin(2, norm=False)
assert(x2[0] == v0)
assert(numpy.sum(x2) == numpy.sum(v))

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
assert(x.npix_allocated == len(a) - 1) # First element was zero
for i in range(0, len(a)):
	assert(x[i] == i)

assert((numpy.asarray(x) == a).all())
assert(x.dense) # Should be dense again

# test ra shifting
x.shift_ra = False
x.ringsparse = True
ki, vi = x.nonzero_pixels()
ii = numpy.argsort(ki)
ki = numpy.asarray(ki)[ii]
vi = numpy.asarray(vi)[ii]

x.shift_ra = True
kr, vr = x.nonzero_pixels()
ii = numpy.argsort(kr)
kr = numpy.asarray(kr)[ii]
vr = numpy.asarray(vr)[ii]
assert((ki == kr).all())
assert((vi == vr).all())

x.shift_ra = False
ki, vi = x.nonzero_pixels()
ii = numpy.argsort(ki)
ki = numpy.asarray(ki)[ii]
vi = numpy.asarray(vi)[ii]
assert((ki == kr).all())
assert((vi == vr).all())


# Conversion to/from flatsky maps
fm_stub = coordinateutils.FlatSkyMap(
        300, 300, core.G3Units.arcmin, proj=coordinateutils.MapProjection.ProjZEA
)
fm = coordinateutils.maputils.healpix_to_flatsky(x, map_stub=fm_stub)
x2 = coordinateutils.maputils.flatsky_to_healpix(fm, map_stub=x.Clone(False))

hitpix = numpy.asarray(x2) > 0
assert(numpy.allclose(numpy.asarray(x)[hitpix], numpy.asarray(x2)[hitpix]))
