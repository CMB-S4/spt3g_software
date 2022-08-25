import numpy as np
from spt3g.core import G3Timestream, G3VectorDouble, DoubleVector, G3VectorInt, IntVector

s = 5.0  # python float
si = 5  # python integer
v = np.abs(np.random.randn(10))  # numpy array
vi = np.random.randint(1, 10, 10)
f = v[0]  # numpy float
fi = vi[0]  # numpy integer

for a in [s, f, v, si, fi, vi, G3Timestream(v), G3VectorDouble(v), DoubleVector(v), G3VectorInt(vi), IntVector(vi)]:
    for b in [s, f, v, si, fi, vi, G3Timestream(v), G3VectorDouble(v), DoubleVector(v)]:
        print(type(a), type(b))
        # binary operators
        v1 = a if np.isscalar(a) or isinstance(a, np.ndarray) else vi if isinstance(a, (G3VectorInt, IntVector)) else v
        v2 = b if np.isscalar(b) or isinstance(b, np.ndarray) else vi if isinstance(b, (G3VectorInt, IntVector)) else v
        assert np.array_equal(a + b, v1 + v2)
        assert np.array_equal(a - b, v1 - v2)
        assert np.array_equal(a * b, v1 * v2)
        assert np.array_equal(a / b, v1 / v2)
        assert np.array_equal(a // b, v1 // v2)
        assert np.array_equal(a % b, v1 % v2)
        assert np.array_equal(a ** b, v1 ** v2)

        if np.isscalar(a) or isinstance(a, (np.ndarray, G3VectorInt, IntVector)):
            continue

        # in-place binary operators
        c = a.__class__(a)
        c += b
        assert isinstance(c, a.__class__)
        assert np.array_equal(c, a + b)

        c = a.__class__(a)
        c -= b
        assert isinstance(c, a.__class__)
        assert np.array_equal(c, a - b)

        c = a.__class__(a)
        c *= b
        assert isinstance(c, a.__class__)
        assert np.array_equal(c, a * b)

        c = a.__class__(a)
        c /= b
        assert isinstance(c, a.__class__)
        assert np.array_equal(c, a / b)

    if isinstance(a, np.ndarray) or np.isscalar(a):
        continue

    for attr in ["sum", "mean", "any", "all", "min", "max", "argmin", "argmax", "var", "std"]:
        assert getattr(a, attr)() == getattr(np, attr)(a)
