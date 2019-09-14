#!/usr/bin/env python
from spt3g import core
import numpy as np

a = core.quat(2,3,4,5)

assert(a+a == 2*a)
assert(a+a == a*2)
assert(a*a == a**2)
assert(a*a*a == a**3)

b = core.G3VectorQuat([a, a**2, 2*a])

assert(b[0].a == 2)
assert(b[0].b == 3)
assert(b[0].c == 4)
assert(b[0].d == 5)

assert(b[1].a == -46)
assert(b[1].b == 12)
assert(b[1].c == 16)
assert(b[1].d == 20)

assert(b[2].c == 8)

c = np.asarray(b)

assert(c.shape == (3,4))

assert(core.quat(*c[0]) == a)
assert(core.quat(*c[1]) == b[1])
assert(core.quat(*c[1]) != b[2])

d = core.G3VectorQuat(c)

assert(d[0] == b[0])
assert(d[1] == b[1])
assert(d[1] != b[2])
assert(d[2] == b[2])

c[1,3] = 15
assert(b[1].d == 15)

e = 2*b
assert((2*np.asarray(b) == np.asarray(e)).all())
e = b/2
assert((np.asarray(b)/2 == np.asarray(e)).all())

assert((np.asarray(b*b) == np.asarray(core.G3VectorQuat([x**2 for x in b]))).all())

