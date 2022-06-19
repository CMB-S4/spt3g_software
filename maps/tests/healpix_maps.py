#!/usr/bin/env python

import numpy as np
from spt3g import core, maps

# Test initialization from sparse arrays

a = np.arange(1500,dtype='float')
a[0] = -1
b = np.arange(1500,dtype='int')
x = maps.HealpixSkyMap((b, a, 64), True, False, maps.MapCoordReference.Equatorial)
x.shift_ra = False
assert(x.nside == 64)
assert(x.npix_allocated == 1500)
assert(x[1499] == 1499)
assert(x[1501] == 0)

# check attributes
attrs = ["nside", "res", "nested", "shift_ra", "units", "pol_type", "pol_conv", "coord_ref", "weighted"]
xc = x.clone()
for attr in attrs:
    assert getattr(xc, attr) == getattr(x, attr), "Attribute {} mismatched on clone()".format(attr)
xc = x.clone(False)
for attr in attrs:
    assert getattr(xc, attr) == getattr(x, attr), "Attribute {} mismatched on clone(False)".format(attr)

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
p0 = hp.nest2ring(64, (hp.ring2nest(32, 0) * 4 + np.arange(4)).astype(int))
v0 = sum([x[int(i)] for i in p0])

x2 = x.rebin(2, norm=False)
assert(x2[0] == v0)
assert(np.sum(x2) == np.sum(v))

# angles
pixels = np.arange(x.size, dtype=int)
alpha, delta = x.pixels_to_angles(pixels)
alpha, delta = np.asarray(alpha), np.asarray(delta)
alpha[alpha < 0] += 360 * core.G3Units.deg
theta, phi = hp.pix2ang(x.nside, pixels)

assert(np.allclose(alpha, phi))
assert(np.allclose(delta, np.pi / 2 - theta))

pixels2 = x.angles_to_pixels(alpha, delta)
assert(np.allclose(pixels, pixels2))

# interpolation
dx = x.res / 2.0
vx = np.asarray(x.get_interp_values(alpha + dx, delta + dx))
vh = hp.get_interp_val(np.asarray(x), theta - dx, phi + dx)
assert(np.allclose(vx, vh))

x.shift_ra = True

# check attributes
xc = x.clone()
for attr in attrs:
    assert getattr(xc, attr) == getattr(x, attr), "Attribute {} mismatched on clone()".format(attr)
xc = x.clone(False)
for attr in attrs:
    assert getattr(xc, attr) == getattr(x, attr), "Attribute {} mismatched on clone(False)".format(attr)

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
assert(np.sum(x2) == np.sum(v))

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
assert(np.sum(x2) == np.sum(v))

x.indexedsparse = True # Dense to indexed
assert(x.nside == 64)
assert(x.npix_allocated == 1500)
for i in range(1,1500):
    assert(x[i] == i)
assert(x[1501] == 0)

# Initialization from dense arrays
a = np.arange(49152)
x = maps.HealpixSkyMap(a)

for i in range(0, len(a)):
    assert(x[i] == i)

assert((np.asarray(x) == a).all())

# Conersion to ring-sparse again (trickiest, this makes sure we get all rings)
x.ringsparse = True
assert(x.npix_allocated == len(a) - 1) # First element was zero
for i in range(0, len(a)):
    assert(x[i] == i)

assert((np.asarray(x) == a).all())
assert(x.dense) # Should be dense again

# test ra shifting
x.shift_ra = False
x.ringsparse = True
ki, vi = x.nonzero_pixels()
ii = np.argsort(ki)
ki = np.asarray(ki)[ii]
vi = np.asarray(vi)[ii]

x.shift_ra = True
kr, vr = x.nonzero_pixels()
ii = np.argsort(kr)
kr = np.asarray(kr)[ii]
vr = np.asarray(vr)[ii]
assert((ki == kr).all())
assert((vi == vr).all())

x.shift_ra = False
ki, vi = x.nonzero_pixels()
ii = np.argsort(ki)
ki = np.asarray(ki)[ii]
vi = np.asarray(vi)[ii]
assert((ki == kr).all())
assert((vi == vr).all())

# Conversion to/from flatsky maps
fm_stub = maps.FlatSkyMap(
    300, 300, core.G3Units.arcmin, proj=maps.MapProjection.ProjZEA
)
fm = maps.maputils.healpix_to_flatsky(x, map_stub=fm_stub)
x2 = maps.maputils.flatsky_to_healpix(fm, map_stub=x.clone(False))

hitpix = np.asarray(x2) > 0
assert(np.allclose(np.asarray(x)[hitpix], np.asarray(x2)[hitpix]))


# Coordinate system rotations
a = np.arange(49152)
x = maps.HealpixSkyMap(a)
alpha, delta = maps.get_ra_dec_map(x)

# equatorial to galactic
xgal = x.clone(False)
x.coord_ref = maps.MapCoordReference.Equatorial
xgal.coord_ref = maps.MapCoordReference.Galactic
maps.reproj_map(x, xgal)
ra, dec = maps.azel.convert_gal_to_radec(alpha, delta)
pix = xgal.angles_to_pixels(ra, dec)
assert(np.allclose(np.asarray(pix), np.asarray(xgal)))

# ... and back
xeq = x.clone(False)
x.coord_ref = maps.MapCoordReference.Galactic
xeq.coord_ref = maps.MapCoordReference.Equatorial
maps.reproj_map(x, xeq)
lon, lat = maps.azel.convert_radec_to_gal(alpha, delta)
pix = xeq.angles_to_pixels(lon, lat)
assert(np.allclose(np.asarray(pix), np.asarray(xeq)))
