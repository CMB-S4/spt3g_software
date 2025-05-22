#!/usr/bin/env python

import numpy
from spt3g import core, maps
from spt3g.maps import FlatSkyMap, MapCoordReference, MapProjection

# 1D arrays
m = FlatSkyMap(1, 500, core.G3Units.arcmin)
assert(m.shape[0] == 500)
m[15] = 65.4
assert(numpy.asarray(m)[15] == 65.4)
assert(numpy.asarray(m).shape == m.shape)

# 2D arrays
m = FlatSkyMap(500, 20, core.G3Units.arcmin)
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
m = FlatSkyMap(w, core.G3Units.arcmin, coord_ref=MapCoordReference.Equatorial)
assert(len(m.shape) == 2)
assert(m.shape[0] == w.shape[0])
assert(m.shape[1] == w.shape[1])
assert(w[17,32] == m[17,32])
assert((w == m).all())

# Test conversion between representations is lossless
# Along the way, test nonzero_pixels() and rebin() once for each
# representation

a = numpy.arange(1500, dtype=float)
a[0] = -1
v = numpy.zeros((50, 50), dtype=float)
v.ravel()[:1500] = a
x = FlatSkyMap(v, core.G3Units.arcmin, proj=MapProjection.ProjZEA)

try:
    v2 = numpy.asarray(x)[10:20, 10:20] # non-contiguous array
    x2 = FlatSkyMap(v2, core.G3Units.arcmin)
    numpy.testing.assert_array_equal(x2, v2)
except Exception as e:
    print(str(e))
    pass
else:
    raise RuntimeError("Expected error!")

assert(not x.sparse)
assert(x.npix_allocated == x.size)
assert(x[1499] == 1499)
assert(x[1501] == 0)

k, v = x.nonzero_pixels()
assert(len(k) == len(v) == 1500)
assert(set(v) == set(a))

v0 = x[0, 0] + x[0, 1] + x[1, 0] + x[1, 1]

x2 = x.rebin(2, norm=False)
assert(x2[0] == v0)
assert(numpy.sum(x2) == numpy.sum(v))
assert(numpy.allclose(x.xy_to_angle(-0.5, -0.5), x2.xy_to_angle(-0.5, -0.5)))

x.sparse = True
assert(x.npix_allocated == 1500)
for i in range(1, 1500):
    assert(x[i] == i)
assert(x[1501] == 0)

k, v = x.nonzero_pixels()
assert(len(k) == len(v) == 1500)
assert(set(v) == set(a))

x2 = x.rebin(2, norm=False)
assert(x2[0] == v0)
assert(numpy.sum(x2) == numpy.sum(v))

x.sparse = False
assert(x.npix_allocated == x.size)
for i in range(1, 1500):
    assert(x[i] == i)
assert(x[1501] == 0)

k, v = x.nonzero_pixels()
assert(len(k) == len(v) == 1500)
assert(set(v) == set(a))

x2 = x.rebin(2, norm=False)
assert(x2[0] == v0)
assert(numpy.sum(x2) == numpy.sum(v))

# query disc
alpha, delta = maps.get_ra_dec_map(x)

from astropy.coordinates import SkyCoord
from astropy import units as u

angles = SkyCoord(
    numpy.asarray(alpha) / core.G3Units.deg * u.degree,
    numpy.asarray(delta) / core.G3Units.deg * u.degree,
).ravel()
radius = 5 * x.res

avec = []
dvec = []
rvec = []
masked = set()

for a in numpy.linspace(numpy.min(alpha), numpy.max(alpha), 20):
    for d in numpy.linspace(numpy.min(delta), numpy.max(delta), 20):
        avec.append(a)
        dvec.append(d)
        rvec.append(radius)
        pix1 = x.query_disc(a, d, radius)
        pix2 = numpy.where(
            angles.separation(
                SkyCoord(
                    a / core.G3Units.deg * u.degree,
                    d / core.G3Units.deg * u.degree,
                )
            ).deg * core.G3Units.deg < radius
        )[0]
        masked |= set(pix2)
        assert not (set(pix1) ^ set(pix2))

mask = maps.make_point_source_mask(x, numpy.asarray(avec), numpy.asarray(dvec), numpy.asarray(rvec))
assert not (set(mask.nonzero()) ^ masked)
