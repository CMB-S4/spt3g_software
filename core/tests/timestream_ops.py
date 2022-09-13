import numpy as np
from spt3g.core import (
    G3Timestream,
    G3VectorDouble,
    DoubleVector,
    G3VectorInt,
    IntVector,
    G3VectorComplexDouble,
    ComplexDoubleVector,
)

s = 5.0  # python float
si = 5  # python integer
sc = s + 1j * s
v = np.abs(np.random.randn(10))  # numpy array
vi = np.random.randint(1, 10, 10)
vc = v + 1j * v  # numpy complex array
f = v[0]  # numpy float
fi = vi[0]  # numpy integer
fc = vc[0]  # numpy complex float

for a in [
    s,
    f,
    v,
    si,
    fi,
    vi,
    sc,
    fc,
    vc,
    G3Timestream(v),
    G3VectorDouble(v),
    DoubleVector(v),
    G3VectorInt(vi),
    IntVector(vi),
    G3VectorComplexDouble(vc),
    ComplexDoubleVector(vc),
]:
    for b in [
        s,
        f,
        v,
        si,
        fi,
        vi,
        sc,
        fc,
        vc,
        G3Timestream(v),
        G3VectorDouble(v),
        DoubleVector(v),
        G3VectorInt(vi),
        IntVector(vi),
        G3VectorComplexDouble(vc),
        ComplexDoubleVector(vc),
    ]:
        print(type(a), type(b))
        # binary operators
        if np.isscalar(a) or isinstance(a, np.ndarray):
            v1 = a
        elif isinstance(a, (G3VectorInt, IntVector)):
            v1 = vi
        elif isinstance(a, (G3VectorComplexDouble, ComplexDoubleVector)):
            v1 = vc
        else:
            v1 = v
        if np.isscalar(b) or isinstance(b, np.ndarray):
            v2 = b
        elif isinstance(b, (G3VectorInt, IntVector)):
            v2 = vi
        elif isinstance(b, (G3VectorComplexDouble, ComplexDoubleVector)):
            v2 = vc
        else:
            v2 = v
        assert np.array_equal(a + b, v1 + v2)
        assert np.array_equal(a - b, v1 - v2)
        assert np.array_equal(a * b, v1 * v2)
        assert np.array_equal(a / b, v1 / v2)
        if not np.iscomplex(b).any() and not np.iscomplex(a).any():
            assert np.array_equal(a // b, v1 // v2)
            assert np.array_equal(a % b, v1 % v2)
        assert np.array_equal(a ** b, v1 ** v2)

        if np.isscalar(a) or isinstance(a, (np.ndarray, G3VectorInt, IntVector)):
            continue

        if np.iscomplex(b).any():
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

    for attr in [
        "sum",
        "mean",
        "any",
        "all",
        "min",
        "max",
        "argmin",
        "argmax",
        "var",
        "std",
    ]:
        assert getattr(a, attr)() == getattr(np, attr)(a)
