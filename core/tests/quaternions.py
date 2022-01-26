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

assert(1./a == core.quat(1.,0,0,0)/a)
assert(np.allclose(core.G3VectorQuat([1./a]), core.G3VectorQuat([(~a) / abs(a)**2])))
assert((np.asarray(abs(b) - abs(~b)) == 0).all())
assert(a/b[0] == (a/b)[0])
assert(b[1]/a == (b/a)[1])

# Test for support of numpy slicing and conversions
quats = np.array([[1., 2., 3., 4., 5.],   # a
                  [0., 0., 0., 3., -1.],   # b
                  [0., 0., 0., 0., 0.],   # c
                  [18., -23., 5., 0., 0.]])  # d

try:
	q = core.G3VectorQuat(quats)
except TypeError:
	# Wrong shape, should fail
	pass
else:
	raise TypeError('Wrong shape, should not convert')

# Non-trivial strides
q = core.G3VectorQuat(quats[:,:4])
assert(q[0] == core.quat(*quats[0,:4]))
assert(q[1] == core.quat(*quats[1,:4]))
assert(q[2] == core.quat(*quats[2,:4]))
assert(q[3] == core.quat(*quats[3,:4]))

# When transposed, has right shape to convert
q = core.G3VectorQuat(quats.T) # Strides, but simple ones
assert(q[0] == core.quat(1,0,0,18))

# Trivial case, packed
q = core.G3VectorQuat(quats.T.copy())
assert(q[0] == core.quat(1,0,0,18))

# Test conversion of integers

qint = np.asarray(quats[:,:4], dtype='int64')
q = core.G3VectorQuat(qint)
assert(q[0] == core.quat(1,2,3,4))
assert(q[1] == core.quat(0,0,0,3))
assert(q[3] == core.quat(18, -23, 5, 0))

