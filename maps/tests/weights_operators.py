#!/usr/bin/env python

import numpy as np
from spt3g import core
from spt3g.maps import G3SkyMapWeights, FlatSkyMap, MapPolType, MapPolConv, MapProjection
from spt3g.maps import remove_weights, remove_weights_t, apply_weights, apply_weights_t, flatten_pol

for pol in [True, False]:
    # allocation
    pol_conv = MapPolConv.IAU if pol else MapPolConv.none
    m = FlatSkyMap(500, 20, core.G3Units.arcmin, pol_conv=pol_conv, proj=MapProjection.Proj5)
    mt = m.clone(False)
    mt.pol_type = MapPolType.T
    if pol:
        mq = m.clone(False)
        mq.pol_type = MapPolType.Q
        mu = m.clone(False)
        mu.pol_type = MapPolType.U
    mw = G3SkyMapWeights(m)
    assert(mw.size == m.size)
    assert(mw.shape == m.shape)
    assert(mw.npix_allocated == 0)
    assert(mw.polarized == pol)
    assert(mw.congruent)

    assert(mw.compatible(m))
    for k in mw.keys():
        assert(mw[k].compatible(mw.TT))

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

    mat = np.array([[2, 0.2, 0.1],
                    [0.2, 0.5, 0.05],
                    [0.1, 0.05, 0.25]]) if pol else 2.
    mw[15] = mat
    assert(mw.npix_allocated == 1)
    assert(np.allclose(mw[15], mat))

    # operators
    mw += 2 * mw
    assert(np.allclose(mw[15], mat * 3))
    assert(mw.npix_allocated == 1)

    mw -= mw / 3
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
    if pol:
        mq *= m
        mu *= m
    mw *= m
    assert(np.allclose(mw[15], mat * 5))
    assert(mw.npix_allocated == 1)

    # weights handling
    if pol:
        remove_weights(mt, mq, mu, mw)
    else:
        remove_weights_t(mt, mw)
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

        iinv = np.linalg.inv(mat * 5)
        inv = mw.inv()
        assert(np.allclose(inv[15], iinv))

    if pol:
        assert(np.isnan([mt[16], mq[16], mu[16]]).all())
    else:
        assert(np.isnan(mt[16]))
    assert(not mt.sparse)
    assert(mt.npix_allocated == mt.size)

    if pol:
        apply_weights(mt, mq, mu, mw)
    else:
        apply_weights_t(mt, mw)
    assert(mt.weighted)
    if pol:
        assert(np.allclose([mt[15], mq[15], mu[15]], vec * 10))
        assert(np.isnan([mt[16], mq[16], mu[16]]).all())
    else:
        assert(np.allclose(mt[15], vec * 10))
        assert(np.isnan(mt[16]))
    assert(not mt.sparse)
    assert(mt.npix_allocated == mt.size)

    # compactify maps back to sparse
    mt.compact(zero_nans=True)
    assert(mt.npix_allocated == 1)
    assert(np.allclose(mt[15], np.atleast_1d(vec * 10)[0]))
    mask = mt.to_mask()
    assert(mask[15] == 1)
    assert(mask.sum() == 1)
    assert(mask.to_map().npix_allocated == 1)

    if pol:
        mq.compact(zero_nans=True)
        mu.compact(zero_nans=True)

    # check memory-efficient weights functions
    if pol:
        remove_weights(mt, mq, mu, mw, zero_nans=True)
    else:
        remove_weights_t(mt, mw, zero_nans=True)
    assert(mt.npix_allocated == 1)

    if pol:
        apply_weights(mt, mq, mu, mw)
        assert(np.allclose([mt[15], mq[15], mu[15]], vec * 10))
    else:
        apply_weights_t(mt, mw)
        assert(np.allclose(mt[15], vec * 10))
    assert(mt.npix_allocated == 1)

    if pol:
        # check flatten_pol
        mq2 = mq.copy()
        mu2 = mu.copy()
        mw2 = mw.copy()
        assert(not mq2.flat_pol)
        assert(not mu2.flat_pol)
        flatten_pol(mq2, mu2, mw2)
        assert(mq2.flat_pol)
        assert(not np.allclose(mq2, mq))
        assert(not np.allclose(mu2, mu))
        assert(not np.allclose(mw2.TQ, mw.TQ))
        assert(not np.allclose(mw2.TU, mw.TU))
        assert(not np.allclose(mw2.QQ, mw.QQ))
        assert(not np.allclose(mw2.QU, mw.QU))
        assert(not np.allclose(mw2.UU, mw.UU))

        flatten_pol(mq2, mu2, mw2, invert=True)
        assert(not mq2.flat_pol)
        assert(np.allclose(mq2, mq))
        assert(np.allclose(mu2, mu))
        assert(np.allclose(mw2.TQ, mw.TQ))
        assert(np.allclose(mw2.TU, mw.TU))
        assert(np.allclose(mw2.QQ, mw.QQ))
        assert(np.allclose(mw2.QU, mw.QU))
        assert(np.allclose(mw2.UU, mw.UU))
