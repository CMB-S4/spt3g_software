import numpy as np
from spt3g.core import G3Timestream, G3VectorDouble, DoubleVector

s = 5.0  # python float
v = np.random.randn(10)  # numpy array
f = v[0]  # numpy float

for a in [s, f, v, G3Timestream(v), G3VectorDouble(v), DoubleVector(v)]:
    for b in [s, f, v, G3Timestream(v), G3VectorDouble(v), DoubleVector(v)]:
        # binary operators
        v1 = a if np.isscalar(a) or isinstance(a, np.ndarray) else v
        v2 = b if np.isscalar(b) or isinstance(b, np.ndarray) else v
        assert np.array_equal(a + b, v1 + v2)
        assert np.array_equal(a - b, v1 - v2)
        assert np.array_equal(a * b, v1 * v2)
        assert np.array_equal(a / b, v1 / v2)

        if np.isscalar(a) or isinstance(a, np.ndarray):
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
        assert getattr(a, attr)() == getattr(v, attr)()
