import numpy as np
from spt3g.core import (
    G3Timestream,
    G3TimestreamMap,
    G3VectorDouble,
    DoubleVector,
    G3VectorInt,
    Int64Vector,
    G3VectorComplexDouble,
    ComplexDoubleVector,
    G3VectorBool,
    BoolVector,
)

s = 5.0  # python float
si = 5  # python integer
sc = s + 1j * s
v = np.abs(np.random.randn(10)) + 1.0  # numpy array
vi = np.random.randint(1, 10, 10)
vc = v + 1j * v  # numpy complex array
f = v[0]  # numpy float
fi = vi[0]  # numpy integer
fc = vc[0]  # numpy complex float
sb = True
vb = v.astype(bool)
fb = vb[0]

all_arr = [
    s,
    f,
    v,
    si,
    fi,
    vi,
    sc,
    fc,
    vc,
    sb,
    fb,
    vb,
    G3Timestream(v),
    G3Timestream(vi),
    G3VectorDouble(v),
    DoubleVector(v),
    G3VectorInt(vi),
    Int64Vector(vi),
    G3VectorComplexDouble(vc),
    ComplexDoubleVector(vc),
    G3VectorBool(vb),
    BoolVector(G3VectorBool(vb)),
]

def check(a1, a2, cls=None, flt=False):
    if flt:
        np.testing.assert_allclose(
            np.asarray(a1),
            np.asarray(a2),
            atol=1e-14,
            rtol=1e-14,
            verbose=True,
        )
    else:
        np.testing.assert_array_equal(np.asarray(a1), np.asarray(a2), verbose=True)
    assert np.asarray(a1).dtype.kind == np.asarray(a2).dtype.kind
    if cls is not None:
        assert isinstance(a1, cls)

def copy(v):
    if isinstance(v, np.ndarray):
        return v.copy()
    return v.__class__(v)

for a in all_arr:
    for b in all_arr:
        n = "/".join([a.__class__.__name__, b.__class__.__name__])
        if not ("Vector" in n or "Timestream" in n):
            continue

        print(type(a), type(b), type(a + b))
        v1 = a if np.isscalar(a) else np.asarray(a)
        v2 = b if np.isscalar(b) else np.asarray(b)
        k1 = np.asarray(a).dtype.kind
        k2 = np.asarray(b).dtype.kind
        flt = k1 in 'fc' or k2 in 'fc'

        # binary operators
        if not (k1 == 'b' and k2 == 'b'):
            check(a + b, v1 + v2, flt=flt)
            check(a - b, v1 - v2, flt=flt)
            check(a * b, v1 * v2, flt=flt)
            check(a / b, v1 / v2, flt=flt)
            if k1 != 'c' and k2 != 'c':
                check(a // b, v1 // v2)
                check(a % b, v1 % v2, flt=flt)
            check(a ** b, v1 ** v2, flt=flt)

        if k1 in 'iub' and k2 in 'iub':
            check(a | b, v1 | v2)
            check(a & b, v1 & v2)
            check(a ^ b, v1 ^ v2)
            check(a << 1, v1 << 1)
            check(a >> 1, v1 >> 1)

        cls = None if np.isscalar(a) or isinstance(a, np.ndarray) else G3VectorBool if 'G3' in n else BoolVector
        check(a < b, v1 < v2, cls)
        check(a <= b, v1 <= v2, cls)
        check(a == b, v1 == v2, cls)
        check(a != b, v1 != v2, cls)
        check(a >= b, v1 >= v2, cls)
        check(a > b, v1 > v2, cls)

        if np.isscalar(a) or k1 == 'b':
            continue

        korder = list('buifc')
        if korder.index(np.asarray(a).dtype.kind) < korder.index(np.asarray(b).dtype.kind):
            continue

        # in-place binary operators
        c = copy(a)
        c += b
        check(c, a + b, a.__class__, flt=flt)

        c = copy(a)
        c -= b
        check(c, a - b, a.__class__, flt=flt)

        c = copy(a)
        c *= b
        check(c, a * b, a.__class__, flt=flt)

        if k1 in 'iub':
            c = copy(a)
            c |= b
            check(c, a | b, a.__class__)

            c = copy(a)
            c &= b
            check(c, a & b, a.__class__)

            c = copy(a)
            c ^= b
            check(c, a ^ b, a.__class__)

            continue

        c = copy(a)
        c /= b
        check(c, a / b, a.__class__, flt=flt)

    if np.isscalar(a) or isinstance(a, np.ndarray):
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
        assert getattr(a, attr)() == getattr(np, attr)(v1)

    if isinstance(a, G3Timestream):
        tsm = G3TimestreamMap(["a", "b", "c"], np.tile(v1, (3, 1)))

        for attr in ["std", "var"]:
            np.testing.assert_allclose(
                getattr(tsm, attr)(axis=-1),
                getattr(np, attr)(np.asarray(tsm), axis=-1),
                atol=1e-14,
                rtol=1e-14,
                verbose=True,
            )
