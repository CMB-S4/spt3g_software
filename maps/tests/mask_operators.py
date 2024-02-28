#!/usr/bin/env python

import numpy as np
from spt3g import core, maps

maplist = [
    maps.FlatSkyMap(np.random.randn(64, 12 * 64), core.G3Units.arcmin),
    maps.HealpixSkyMap(np.random.randn(12 * 64 * 64)),
]

for x in maplist:
    print(x)
    m = x.to_mask()

    # attributes
    assert(m.size == x.size)
    assert(m.shape == x.shape)

    # element-wise assignment
    if isinstance(x, maps.FlatSkyMap):
        v = m[5, 5]
        m[5, 5] = not v
        assert(m[5, 5] != v)
        m[5, 5] = v
        assert(m[5, 5] == v)
    else:
        v = m[5]
        m[5] = not v
        assert(m[5] != v)
        m[5] = v
        assert(m[5] == v)

    # masked assignment and retrieval
    print("assignment")
    pospixels = x[x > 0]
    assert(len(pospixels) == (x > 0).sum())
    assert((np.asarray(pospixels) > 0).all())

    x2 = x.__class__(x)
    x2[x2 > 0] = 0
    assert((x2 > 0).sum() == 0)
    nzero = (x2 == 0).sum()
    x2[x2 == 0] = np.abs(np.random.randn((x2 == 0).sum())) + 1
    assert((x2 > 0).sum() == nzero)
    assert((x2 == 0).sum() == 0)

    # init and clone
    mx = maps.G3SkyMapMask(x, True)
    assert((mx == m).all())
    mn = maps.G3SkyMapMask(x, False)
    assert(mn.sum() == 0)
    a = np.asarray(x).ravel()
    ma = m.array_clone(a)
    assert((m == ma).all())
    mp = m.array_clone(a > 0)
    assert(((x > 0) == mp).all())

    # float comparison
    print("float comp")
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
    print("map comp")
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
    print("logic")
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
    print("conversion")
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

    mpix = m1.nonzero()
    xpix, _ = x4.nonzero_pixels()
    xpixm = x4.nonzero()

    assert(len(set(mpix) ^ set(xpix)) == 0)
    assert(len(set(mpix) ^ set(xpixm)) == 0)

    m1.apply_mask(m2)
    assert(not m1.any())

    # Direct array conversion
    assert(np.asarray(m2).shape == np.asarray(m2.to_map()).shape)
    assert((np.asarray(m2) == np.asarray(m2.to_map())).all())

    # ufuncs
    print("ufuncs")
    if isinstance(x, maps.FlatSkyMap):
        x.sparse = True
    else:
        x.ringsparse = True
    m = x < 0

    for attr in ["all", "any", "sum", "mean", "var", "std", "min", "max", "argmin", "argmax"]:
        for w in [None, m]:
            if w is not None and attr in ["argmin", "argmax"]:
                continue
            kwargs = {"where": w} if w is not None else {}
            if attr in ["std", "var"]:
                kwargs["ddof"] = 1
            try:
                v1 = getattr(np, attr)(x, **kwargs)
            except TypeError as e:
                # ignore errors with older numpy versions
                if w is not None and "unexpected keyword argument 'where'" in str(e):
                    continue
                raise
            v2 = getattr(x, attr)(**kwargs)
            if w is not None:
                kwargs["where"] = np.asarray(w.to_map(), dtype=bool)
                if attr in ["min", "max"]:
                    kwargs["initial"] = np.inf if attr == "min" else -np.inf
            v3 = getattr(np, attr)(np.asarray(x.copy()), **kwargs)
            print(attr, w, v1, v2, v3)
            assert(v1 == v2)
            assert(np.isclose(v2, v3))
            if isinstance(x, maps.FlatSkyMap):
                assert(x.sparse == True)
            else:
                assert(x.ringsparse == True)

    x[m.nonzero()[0]] = np.nan

    for attr in ["nansum", "nanmean", "nanvar", "nanstd", "nanmin", "nanmax", "nanargmin", "nanargmax"]:
        v1 = getattr(np, attr)(x)
        v2 = getattr(x, attr)()
        v3 = getattr(np, attr)(np.asarray(x.copy()))
        print(attr, v1, v2, v3)
        assert(np.isclose(v1, v2))
        assert(np.isclose(v2, v3))

    x[m.nonzero()[1]] = np.inf

    for attr in ["isinf", "isnan", "isfinite"]:
        for w in [None, m]:
            print(attr, w)
            kwargs = {"where": w} if w is not None else {}
            m1 = getattr(x, attr)(**kwargs)
            if w is not None:
                kwargs["where"] = np.asarray(w.to_map(), dtype=bool)
                kwargs["out"] = np.zeros_like(kwargs["where"])
            m2 = maps.G3SkyMapMask(x, getattr(np, attr)(x, **kwargs).ravel())
            assert (m1 == m2).all()
