#!/usr/bin/env python

import numpy as np
from spt3g import core, maps

shape = (300, 500)
x = maps.FlatSkyMap(np.random.randn(*shape), core.G3Units.arcmin)
m = x.to_mask()

# attributes
assert(m.size == x.size)
assert(m.shape == x.shape)

# element-wise assignment
v = m[5, 5]
m[5, 5] = not v
assert(m[5, 5] != v)
m[5, 5] = v
assert(m[5, 5] == v)

# float comparison
m1 = x == 0
m2 = x != 0
assert(m1.sum() + m2.sum() == x.size)

m1 = x <= 0
m2 = x > 0
assert(m1.sum() + m2.sum() == x.size)

m1 = x < 0
m2 = x >= 0
assert(m1.sum() + m2.sum() == x.size)

# map comparison
y = 2 * x.clone()

m1 = y == x
m2 = y != x
assert(m1.sum() + m2.sum() == x.size)

m1 = y <= x
m2 = y > x
assert(m1.sum() + m2.sum() == x.size)

m1 = y < x
m2 = y >= x
assert(m1.sum() + m2.sum() == x.size)

# logic operators
assert((m1 == ~m2).all())
assert(not (m1 & m2).any())
assert((m1 | m2).all())
assert((m1 ^ m2).all())

m3 = m1.clone()
m3 |= m2
assert(m3.all())
m3 = m1.clone()
m3 &= m2
assert(not m3.any())
m3 = m1.clone()
m3 ^= m2
assert(m3.all())
m3 = m1.clone()
m3.invert()
assert((m3 == m2).all())

# convert to mask
x2 = x.clone()
x2 *= m1
assert((x2.to_mask() == m1).all())
assert(((x * m1).to_mask() == m1).all())

x3 = x.clone()
x3.apply_mask(m1, inverse=True)
assert((x3.to_mask() == m2).all())

x4 = m1.to_map()
assert(np.sum(x4) == m1.sum())
assert((x4.to_mask() == m1).all())

mpix = m1.nonzero_pixels()
xpix, _ = x4.nonzero_pixels()

assert(len(set(mpix) ^ set(xpix)) == 0)

m1.apply_mask(m2)
assert(not m1.any())
