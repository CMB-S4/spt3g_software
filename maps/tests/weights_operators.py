#!/usr/bin/env python

import numpy as np
from spt3g import core
from spt3g.maps import G3SkyMapWeights, FlatSkyMap, MapPolType, MapPolConv
from spt3g.maps import get_mask_map, remove_weights, apply_weights

for pol in [True, False]:
    # allocation
    m = FlatSkyMap(500, 20, core.G3Units.arcmin, pol_conv=MapPolConv.IAU)
    mt = m.Clone(False)
    mt.pol_type = MapPolType.T
    if pol:
        mq = m.Clone(False)
        mq.pol_type = MapPolType.Q
        mu = m.Clone(False)
        mu.pol_type = MapPolType.U
    mw = G3SkyMapWeights(m, polarized=pol)
    assert(mw.size == m.size)
    assert(mw.shape == m.shape)
    assert(mw.npix_allocated == 0)
    assert(mw.polarized == pol)
    assert(mw.congruent)

    assert(mw.IsCompatible(m))
    assert(mw['TT'].IsCompatible(mw.TT))

    if not pol:
        assert(mw.QQ is None)
        assert(mw.UU is None)

    # assignment
    vec = np.array([10., 1., 1.]) if pol else 10.
    if pol:
        mt[15], mq[15], mu[15] = vec
        assert(np.allclose([mt[15], mq[15], mu[15]], vec))
    else:
        mt[15] = vec
        assert(np.allclose(mt[15], vec))

    mat = np.array([[2, 0.1, 0.1],
                    [0.1, 0.5, 0.05],
                    [0.1, 0.05, 0.5]]) if pol else 2.
    mw[15] = mat
    assert(mw.npix_allocated == 1)
    assert(np.allclose(mw[15], mat))

    # operators
    mw += mw
    assert(np.allclose(mw[15], mat * 2))
    assert(mw.npix_allocated == 1)

    mw *= 2
    assert(np.allclose(mw[15], mat * 4))
    assert(mw.npix_allocated == 1)

    mw /= 4
    assert(np.allclose(mw[15], mat))

    mw /= 2
    assert(np.allclose(mw[15], mat / 2))

    m[15] = 10
    mt *= m
    mq *= m
    mu *= m
    mw *= m
    assert(np.allclose(mw[15], mat * 5))
    assert(mw.npix_allocated == 1)

    # weights handling

    remove_weights(mt, mq, mu, mw)
    assert(not mt.weighted)
    assert(np.allclose(mw[15], mat * 5))
    ivec = np.linalg.solve(mat * 5, vec * 10) if pol else vec * 2 / mat
    if pol:
        assert(np.allclose([mt[15], mq[15], mu[15]], ivec))
    else:
        assert(np.allclose(mt[15], ivec))

    if pol:
        idet = np.linalg.det(mat * 5)
        det = mw.det()
        assert(np.allclose(det[15], idet))

        icond = np.linalg.cond(mat * 5)
        cond = mw.cond()
        assert(np.allclose(cond[15], icond))

    if pol:
        assert(np.isnan([mt[16], mq[16], mu[16]]).all())
    else:
        assert(np.isnan(mt[16]))
    assert(not mt.sparse)
    assert(mt.npix_allocated == mt.size)

    apply_weights(mt, mq, mu, mw)
    assert(mt.weighted)
    if pol:
        assert(np.allclose([mt[15], mq[15], mu[15]], vec * 10))
        assert(np.isnan([mt[16], mq[16], mu[16]]).all())
    else:
        assert(np.allclose(mt[15], vec * 10))
        assert(np.isnan(mt[16]))
    assert(not mt.sparse)
    assert(mt.npix_allocated == mt.size)

    mt.compact(zero_nans=True)
    assert(mt.npix_allocated == 1)
    assert(np.allclose(mt[15], np.atleast_1d(vec * 10)[0]))
    mask = get_mask_map(mt)
    assert(mask[15] == 1)
    assert(mask.npix_allocated == 1)

    # compactify maps back to sparse
    if pol:
        mq.compact(zero_nans=True)
        mu.compact(zero_nans=True)

    # check memory-efficient weights functions
    remove_weights(mt, mq, mu, mw, zero_nans=True)
    assert(mt.npix_allocated == 1)

    apply_weights(mt, mq, mu, mw)
    assert(mt.npix_allocated == 1)
    assert(mt[15] == np.atleast_1d(vec * 10)[0])
