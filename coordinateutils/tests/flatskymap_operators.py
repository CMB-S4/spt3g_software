#!/usr/bin/env python

import numpy
from spt3g import core
from spt3g.coordinateutils import FlatSkyMap, MapCoordReference

# Simple in-place operators
m = FlatSkyMap(500, 20, core.G3Units.arcmin)
m[15] = 10
assert(m.sparse)

m *= 5
m /= 25
assert(m[15] == 2)
assert(m.sparse)
assert(m[16] == 0)

assert((-m).sparse)
assert((-m)[15] == -2)

m += 3
m -= 14
assert(m[15] == -9)
assert(not m.sparse)
assert(m[16] == -11)

a = -11*numpy.ones(m.shape)
a[0,15] += 2
assert((m == a).all()) # Implicitly tests numpy conversions too

assert((m*0).npix_allocated == 0)

# Map-by-map operations, with two sparse maps, one dense and one sparse,
# and two dense

m += 11
m.sparse = True
assert(m.npix_allocated == 1)

m *= 2 # Get numbers bigger

m1 = m
m2 = m.Clone(True)
m2.sparse = False

for pair in [(m1, m1), (m1, m2), (m2, m1), (m2, m2)]:
	t = pair[0]*pair[1]
	assert(t.sparse == pair[0].sparse)
	assert(t.npix_allocated == pair[0].npix_allocated)
	assert(t[15] == 16)

	t = pair[0]+pair[1]
	assert(t.sparse == pair[0].sparse)
	assert(t.npix_allocated == pair[0].npix_allocated)
	assert(t[15] == 8)

	t = pair[0]-2*pair[1]
	assert(t.sparse == pair[0].sparse)
	assert(t.npix_allocated == pair[0].npix_allocated)
	assert(t[15] == -4)

	t = pair[0]/pair[1]
	assert(t.sparse == pair[0].sparse)
	assert(t.npix_allocated == t.size)
	assert(t[15] == 1)
	assert(not numpy.isfinite(t[12]))

# With a null map
m3 = m.Clone(False)

for pair in [(m1, m3), (m2, m3), (m3, m2), (m3, m1)]:
	nonnull = pair[1] if pair[0] is m3 else pair[0]

	t = pair[0]*pair[1]
	assert(t.sparse == True)
	assert(t.npix_allocated == 0)

	t = pair[0]+pair[1]
	assert(t.sparse == nonnull.sparse)
	assert(t.npix_allocated == nonnull.npix_allocated)
	assert(t[15] == 4)

	t = pair[0]-pair[1]
	assert(t.sparse == nonnull.sparse)
	assert(t.npix_allocated == nonnull.npix_allocated)
	assert(t[15] == -4 or t[15] == 4)

	t = pair[0]/pair[1]
	assert(t.npix_allocated == t.size)
	assert(not numpy.isfinite(t[12]))
	if pair[0] is m3:
		assert(t[15] == 0)
	

