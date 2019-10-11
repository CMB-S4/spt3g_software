#!/usr/bin/env python

import numpy as np
from spt3g import core
from spt3g.coordinateutils import HealpixSkyMap

# Simple in-place operators
m = HealpixSkyMap(64)
m[15] = 10
assert(m.ringsparse)

m *= 5
m /= 25
assert(m[15] == 2)
assert(m.ringsparse)
assert(m[16] == 0)

assert((-m).ringsparse)
assert((-m)[15] == -2)

m += 3
m -=14
assert(m[15] == -9)
assert(not m.ringsparse)
assert(m[16] == -11)

n = 1 - m
assert(n[15] == 10)
assert(n[16] == 12)

a = -11 * np.ones(m.shape)
a[15] += 2
assert((m == a).all()) # Implicitly tests numpy conversions too

assert((m * 0).npix_allocated == 0)

m += 11
m.indexedsparse = True
assert(m.npix_allocated == 1)

n = 2. / m
assert(n[15] == 1)
assert(np.isinf(n[16]))
assert(n.npix_allocated == n.size)

# Map-by-map operations, with pairs of maps of any kind of density

m *= 2 # Get numbers bigger

m1 = m
m2 = m.Clone(True)
m2.ringsparse = True
m3 = m.Clone(True)
m3.dense = True

nm1  = m1.npix_allocated
nm2 = m2.npix_allocated
nm3 = m3.npix_allocated

def check_state(a, b, check_size=True):
    assert(m1.npix_allocated == nm1)
    assert(m2.npix_allocated == nm2)
    assert(m3.npix_allocated == nm3)
    assert(a.indexedsparse == b.indexedsparse)
    assert(a.ringsparse == b.ringsparse)
    assert(a.dense == b.dense)
    if check_size:
        assert(a.npix_allocated == b.npix_allocated)

for pair in [(m1, m1), (m1, m2), (m1, m3),
             (m2, m1), (m2, m2), (m2, m3),
             (m3, m1), (m3, m2), (m3, m3)]:
    t = pair[0] * pair[1]
    check_state(pair[0], t)
    assert(t[15] == 16)

    t = pair[0] + pair[1]
    check_state(pair[0], t)
    assert(t[15] == 8)

    t = pair[0] - 2 * pair[1]
    check_state(pair[0], t)
    assert(t[15] == -4)

    t = pair[0] / pair[1]
    check_state(pair[0], t, False)
    assert(t.npix_allocated == t.size)
    assert(t[15] == 1)
    assert(not np.isfinite(t[12]))

# With a null map
m4 = m.Clone(False)

for pair in [(m1, m4), (m2, m4), (m3, m4),
             (m4, m1), (m4, m2), (m4, m3)]:
    nonnull = pair[1] if pair[0] is m4 else pair[0]

    t = pair[0] * pair[1]
    assert(t.dense == False)
    assert(t.npix_allocated == 0)
    assert(m4.npix_allocated == 0)

    t = pair[0] + pair[1]
    check_state(t, nonnull)
    assert(t[15] == 4)
    assert(m4.npix_allocated == 0)

    t = pair[0] - pair[1]
    check_state(t, nonnull)
    assert(t[15] in [-4, 4])
    assert(m4.npix_allocated == 0)

    t = pair[0] / pair[1]
    assert(t.npix_allocated == t.size)
    assert(not np.isfinite(t[12]))
    if pair[0] is m4:
        assert(t[15] == 0)
    assert(m4.npix_allocated == 0)
