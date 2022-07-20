#!/usr/bin/env python

import numpy as np
from spt3g import core
from spt3g.maps import FlatSkyMap, MapProjection, get_ra_dec_map, get_map_moments, get_map_median
from scipy.stats import skew, kurtosis

# Sparse extension operators
m = FlatSkyMap(500, 20, core.G3Units.arcmin)
m[7345] = 4
m[7345-500] = 4 # Backward
m[7345+500] = 4 # Forward
assert(m.sparse)
assert(m.npix_allocated == 3)
assert(m.npix_nonzero == 3)
m[7345-4*500] = 4 # Several steps back
assert(m.npix_allocated == 6)
assert(m.npix_nonzero == 4)
m[7345+3*500] = 4 # Several steps forward
assert(m.npix_allocated == 8)
assert(m.npix_nonzero == 5)

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

n = 1 - m
assert(n[15] == 10)
assert(n[16] == 12)

a = -11 * np.ones(m.shape)
a[0, 15] += 2
assert((m == a).all()) # Implicitly tests numpy conversions too

assert((m*0).npix_allocated == 0)

m += 11
m.sparse = True
assert(m.npix_allocated == 1)

n = 2. / m
assert(n[15] == 1)
assert(np.isinf(n[16]))
assert(n.npix_allocated == n.size)

# compactification
np.asarray(n)[np.isinf(n)] = np.nan
n.compact(zero_nans=True)
assert(n[16] == 0)
assert(n.npix_allocated == 1)

np.asarray(n)[np.isinf(n)] = np.nan
n.sparse = True
n.compact(zero_nans=True)
assert(n[16] == 0)
assert(n.npix_allocated == 1)

n = m ** 2.
assert(n[15] == 4)
assert(n.npix_allocated == 1)

# Map-by-map operations, with two sparse maps, one dense and one sparse,
# and two dense

m *= 2 # Get numbers bigger
assert((m == n).all())
assert((m > 0).any())
assert((m > 0).sum() == 1)
assert((m > 0).to_map().npix_allocated == 1)

m1 = m
m2 = m.copy()
m2.sparse = False

nm1 = m1.npix_allocated
nm2 = m2.npix_allocated

for pair in [(m1, m1), (m1, m2), (m2, m1), (m2, m2)]:
    t = pair[0] * pair[1]
    assert(t.sparse == pair[0].sparse)
    assert(t.npix_allocated == pair[0].npix_allocated)
    assert(t[15] == 16)
    assert(m1.npix_allocated == nm1)
    assert(m2.npix_allocated == nm2)

    t = pair[0] + pair[1]
    assert(t.sparse == pair[0].sparse)
    assert(t.npix_allocated == pair[0].npix_allocated)
    assert(t[15] == 8)
    assert(m1.npix_allocated == nm1)
    assert(m2.npix_allocated == nm2)

    t = pair[0] - 2 * pair[1]
    assert(t.sparse == pair[0].sparse)
    assert(t.npix_allocated == pair[0].npix_allocated)
    assert(t[15] == -4)
    assert(m1.npix_allocated == nm1)
    assert(m2.npix_allocated == nm2)

    t = pair[0] / pair[1]
    assert(t.sparse == pair[0].sparse)
    assert(t.npix_allocated == t.size)
    assert(t[15] == 1)
    assert(not np.isfinite(t[12]))
    assert(m1.npix_allocated == nm1)
    assert(m2.npix_allocated == nm2)

# With a null map
m3 = m.clone(False)

for pair in [(m1, m3), (m2, m3), (m3, m2), (m3, m1)]:
    nonnull = pair[1] if pair[0] is m3 else pair[0]

    t = pair[0] * pair[1]
    assert(t.sparse == True)
    assert(t.npix_allocated == 0)
    assert(m1.npix_allocated == nm1)
    assert(m2.npix_allocated == nm2)
    assert(m3.npix_allocated == 0)

    t = pair[0] + pair[1]
    assert(t.sparse == nonnull.sparse)
    assert(t.npix_allocated == nonnull.npix_allocated)
    assert(t[15] == 4)
    assert(m1.npix_allocated == nm1)
    assert(m2.npix_allocated == nm2)
    assert(m3.npix_allocated == 0)

    t = pair[0] - pair[1]
    assert(t.sparse == nonnull.sparse)
    assert(t.npix_allocated == nonnull.npix_allocated)
    assert(t[15] == -4 or t[15] == 4)
    assert(m1.npix_allocated == nm1)
    assert(m2.npix_allocated == nm2)
    assert(m3.npix_allocated == 0)

    t = pair[0] / pair[1]
    assert(not np.isfinite(t[12]))
    if pair[0] is m3:
        assert(t[15] == 0)
    assert(m1.npix_allocated == nm1)
    assert(m2.npix_allocated == nm2)
    assert(m3.npix_allocated == 0)


for shape in [(20, 500), (21, 501)]:
    print('shape', shape)

    # patch extraction / insertion
    m = FlatSkyMap(shape[1], shape[0], core.G3Units.arcmin, proj=MapProjection.ProjZEA)
    np.asarray(m)[:] = np.random.randn(*m.shape)
    malpha, mdelta = get_ra_dec_map(m)

    x0 = 45
    y0 = 13

    for dy, dx in [[10, 50], [11, 51]]:
        print('    patch', (dy, dx))
        p = m.extract_patch(45, 13, dx, dy)
        palpha, pdelta = get_ra_dec_map(p)
        x1 = x0 - dx // 2
        y1 = y0 - dy // 2
        sx = slice(x1, x1 + dx)
        sy = slice(y1, y1 + dy)
        assert(np.allclose(np.asarray(malpha)[sy, sx], palpha))
        assert(np.allclose(np.asarray(mdelta)[sy, sx], pdelta))
        assert(np.allclose(np.asarray(m)[sy, sx], p))

        m2 = m.clone(False)
        m2.insert_patch(p)
        assert(np.allclose(np.asarray(m)[sy, sx], np.asarray(m2)[sy, sx]))

    # Slice operators: make sure they work like numpy slicing
    assert((np.asarray(m[10:17,320:482]) == np.asarray(m.copy())[10:17,320:482]).all())
    # But give the right type...
    assert(m[10:17,320:482].__class__ == m.__class__)

    # Try setting things
    old_chunk = m[10:17,320:482]
    m[10:17,320:482] = old_chunk*2
    assert((np.asarray(m.copy())[10:17,320:482] == np.asarray(old_chunk)*2).all())
    m[10:17,320:482] = np.asarray(old_chunk*3)
    assert((np.asarray(m.copy())[10:17,320:482] == np.asarray(old_chunk)*3).all())
    m[10:17,320:482] = old_chunk
    mcopy = m.copy()
    m[:] = np.asarray(m * 2)
    assert((np.asarray(m) == np.asarray(mcopy * 2)).all())
    m[:] = np.asarray(mcopy)

    # Make sure inserting it in the wrong place (where coordinates don't make sense, but numpy
    # would allow it) fails
    failed = False
    try:
        m[11:18,320:482] = old_chunk
    except ValueError:
        failed = True
    assert(failed)

    # negative slice indices
    assert((np.asarray(m.copy())[-10:-3, -180:-18] == np.asarray(m.copy())[-10:-3, -180:-18]).all())

    # padding / cropping, with even and odd changes in dimension
    pad = 10
    for off in [0, 1]:
        print('    padding', 2 * pad + off)

        mpad = m.reshape(m.shape[1] + 2 * pad + off, m.shape[0] + 2 * pad + off)
        assert(mpad.npix_allocated == m.npix_allocated)
        a0 = np.array([m.alpha_center, m.delta_center])
        a1 = np.array([mpad.alpha_center, mpad.delta_center])
        assert(np.allclose(a0, a1))
        x0 = np.array([m.y_center, m.x_center])
        x1 = np.array([mpad.y_center, mpad.x_center])
        v0 = m[m.angle_to_pixel(*a0)]
        v1 = mpad[mpad.angle_to_pixel(*a1)]
        assert(v0 == v1)
        palpha, pdelta = get_ra_dec_map(mpad)
        sx = slice(mpad.shape[1] // 2 - m.shape[1] // 2, mpad.shape[1] // 2 - m.shape[1] // 2 + m.shape[1])
        sy = slice(mpad.shape[0] // 2 - m.shape[0] // 2, mpad.shape[0] // 2 - m.shape[0] // 2 + m.shape[0])
        assert(np.allclose(np.asarray(palpha)[sy, sx], malpha))
        assert(np.allclose(np.asarray(pdelta)[sy, sx], mdelta))
        assert(np.allclose(np.asarray(mpad)[sy, sx], np.asarray(m)))

        mcrop = mpad.reshape(m.shape[1], m.shape[0])
        calpha, cdelta = get_ra_dec_map(mcrop)
        a2 = np.array([mcrop.alpha_center, mcrop.delta_center])
        assert(np.allclose(a0, a2))
        x2 = np.array([mcrop.y_center, mcrop.x_center])
        assert(np.allclose(x0, x2))
        v2 = mcrop[mcrop.angle_to_pixel(*a2)]
        assert(v0 == v2)
        assert(np.allclose(calpha, malpha))
        assert(np.allclose(cdelta, mdelta))
        assert(np.allclose(mcrop, m))
        assert(m.compatible(mcrop))

        # statistics
        m1 = np.asarray(m).ravel()
        stats0 = [np.mean(m1), np.var(m1), skew(m1), kurtosis(m1)]
        stats1 = get_map_moments(m, order=4)
        assert(np.allclose(stats1, stats0))
        med0 = np.median(m1)
        med1 = get_map_median(m)
        assert(np.allclose(med1, med0))

        stats2 = get_map_moments(mpad, order=4, ignore_zeros=True)
        assert(np.allclose(stats2, stats0))
        med2 = get_map_median(mpad, ignore_zeros=True)
        assert(np.allclose(med2, med0))

        np.asarray(mpad)[np.asarray(mpad) == 0] = np.nan
        stats3 = get_map_moments(mpad, order=4, ignore_nans=True)
        assert(np.allclose(stats3, stats0))
        med3 = get_map_median(mpad, ignore_nans=True)
        assert(np.allclose(med3, med0))

# convolution
from scipy.signal import convolve2d
from spt3g.maps import convolve_map

kernel = FlatSkyMap(np.random.randn(5, 5), m.res)
mconv = convolve2d(m, kernel, mode='same')
mconv2 = convolve_map(m, kernel)
assert(np.allclose(mconv, mconv2))
