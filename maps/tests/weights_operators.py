#!/usr/bin/env python

import numpy as np
from spt3g import core
from spt3g.maps import G3SkyMapWithWeights, G3SkyMapWeights, FlatSkyMap
from spt3g.maps import get_mask_map, remove_weights, apply_weights

for pol in [True, False]:
    # allocation
    m = FlatSkyMap(500, 20, core.G3Units.arcmin)
    mw = G3SkyMapWithWeights(m, weighted=True, polarized=pol)
    assert(mw.size == m.size)
    assert(mw.weights.size == m.size)
    assert(mw.shape == m.shape)
    assert(mw.weights.shape == m.shape)
    assert(mw.npix_allocated == 0)
    assert(mw.polarized == pol)
    assert(mw.weights.polarized == pol)
    assert(mw.weighted)
    assert(mw.congruent)

    assert(np.all(mw['T'].IsCompatible(mw.T)))
    assert(np.all(mw.weights['TT'].IsCompatible(mw.weights.TT)))

    if not pol:
        assert(mw.Q is None)
        assert(mw.U is None)
        assert(mw.weights.QQ is None)
        assert(mw.weights.UU is None)

    # assignment
    vec = np.array([10., 1., 1.]) if pol else 10.
    mw[15] = vec
    assert(mw.npix_allocated == 1)
    assert(np.allclose(mw[15], vec))

    mat = np.array([[2, 0.1, 0.1],
                    [0.1, 0.5, 0.05],
                    [0.1, 0.05, 0.5]]) if pol else 2.
    mw.weights[15] = mat
    assert(np.allclose(mw.weights[15], mat))

    # operators
    mw += mw
    assert(np.allclose(mw[15], vec * 2))
    assert(np.allclose(mw.weights[15], mat * 2))
    assert(mw.npix_allocated == 1)

    mw *= 2
    assert(np.allclose(mw[15], vec * 4))
    assert(np.allclose(mw.weights[15], mat * 4))
    assert(mw.npix_allocated == 1)

    mw /= 4
    assert(np.allclose(mw[15], vec))
    assert(np.allclose(mw.weights[15], mat))

    mw.weights /= 2
    assert(np.allclose(mw.weights[15], mat / 2))

    m[15] = 10
    mw *= m
    assert(np.allclose(mw[15], vec * 10))
    assert(np.allclose(mw.weights[15], mat * 5))
    assert(mw.npix_allocated == 1)

    # weights handling
    weights = mw.remove_weights()
    assert(not mw.weighted)
    assert(np.allclose(weights[15], mat * 5))
    ivec = np.linalg.solve(mat * 5, vec * 10) if pol else vec * 2 / mat
    assert(np.allclose(mw[15], ivec))

    if pol:
        idet = np.linalg.det(mat * 5)
        det = weights.det()
        assert(np.allclose(det[15], idet))

        icond = np.linalg.cond(mat * 5)
        cond = weights.cond()
        assert(np.allclose(cond[15], icond))

    assert(np.isnan(mw[16]).all())
    assert(not mw.sparse)
    assert(mw.npix_allocated == mw.size)

    mw.apply_weights(weights)
    assert(mw.weighted)
    assert(np.allclose(mw[15], vec * 10))
    assert(np.isnan(mw[16]).all())
    assert(not mw.sparse)
    assert(mw.npix_allocated == mw.size)

    tmap = mw.T.copy()
    tmap.compact(zero_nans=True)
    assert(tmap.npix_allocated == 1)
    assert(np.allclose(tmap[15], np.atleast_1d(vec * 10)[0]))
    mask = get_mask_map(tmap)
    assert(mask[15] == 1)
    assert(mask.npix_allocated == 1)

    # compactify maps back to sparse
    if pol:
        qmap = mw.Q.copy()
        qmap.compact(zero_nans=True)
        umap = mw.U.copy()
        umap.compact(zero_nans=True)
    else:
        qmap = umap = None

    # check memory-efficient weights functions
    remove_weights(tmap, qmap, umap, weights, zero_nans=True)
    assert(tmap.npix_allocated == 1)

    apply_weights(tmap, qmap, umap, weights)
    assert(tmap.npix_allocated == 1)
    assert(tmap[15] == mw.T[15])
