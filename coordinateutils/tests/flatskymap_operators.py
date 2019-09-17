#!/usr/bin/env python

import numpy as np
from spt3g import core
from spt3g.coordinateutils import FlatSkyMap

# Sparse extension operators
m = FlatSkyMap(500, 20, core.G3Units.arcmin)
m[7345] = 4
m[7345-500] = 4 # Backward
m[7345+500] = 4 # Forward
assert(m.sparse)
assert(m.npix_allocated == 3)
m[7345-4*500] = 4 # Several steps back
assert(m.npix_allocated == 6)
m[7345+3*500] = 4 # Several steps forward
assert(m.npix_allocated == 8)

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

# Map-by-map operations, with two sparse maps, one dense and one sparse,
# and two dense

m *= 2 # Get numbers bigger

m1 = m
m2 = m.Clone(True)
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
m3 = m.Clone(False)

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
    assert(t.npix_allocated == t.size)
    assert(not np.isfinite(t[12]))
    if pair[0] is m3:
        assert(t[15] == 0)
    assert(m1.npix_allocated == nm1)
    assert(m2.npix_allocated == nm2)
    assert(m3.npix_allocated == 0)
