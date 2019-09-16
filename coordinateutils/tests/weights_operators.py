#!/usr/bin/env python

import numpy as np
from spt3g import core
from spt3g.coordinateutils import G3SkyMapWithWeights, G3SkyMapWeights, FlatSkyMap

# allocation
m = FlatSkyMap(500, 20, core.G3Units.arcmin)
mw = G3SkyMapWithWeights(m, isweighted=True, ispolarized=True)
assert(mw.shape == m.shape)
assert(mw.npix_allocated == 0)
assert(mw.polarized)
assert(mw.weighted)

# assignment
vec = np.array([10., 1., 1.])
mw[15] = vec
assert(mw.npix_allocated == 1)
assert(np.allclose(mw[15], vec))

mat = np.array([[2, 0.1, 0.1],
                [0.1, 0.5, 0.05],
                [0.1, 0.05, 0.5]])
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
ivec = np.linalg.solve(mat * 5, vec * 10)
assert(np.allclose(mw[15], ivec))
assert(mw.npix_allocated == 1)

mw.apply_weights(weights)
assert(mw.weighted)
assert(np.allclose(mw[15], vec * 10))
assert(mw.npix_allocated == 1)
