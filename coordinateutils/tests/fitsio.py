#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
from spt3g import core, coordinateutils

arr = np.arange(300 * 300).reshape(300, 300)

for p in [0, 1, 2, 4, 5, 6, 7, 9]:
    print('Checking Proj{}'.format(p))

    proj = getattr(coordinateutils.MapProjection, 'Proj{}'.format(p))
    fm1 = coordinateutils.FlatSkyMap(arr, core.G3Units.arcmin, proj=proj)

    coordinateutils.maputils.save_skymap_fits('test_map.fits', fm1, overwrite=True)
    fm2 = coordinateutils.maputils.load_skymap_fits('test_map.fits')['T']

    assert(fm1.IsCompatible(fm2))
    assert(np.allclose(fm1, fm2))

print('Checking Healpix')
hm1 = coordinateutils.HealpixSkyMap(np.arange(12 * 64 * 64))
coordinateutils.maputils.save_skymap_fits('test_map.fits', hm1, overwrite=True)
hm2 = coordinateutils.maputils.load_skymap_fits('test_map.fits')['T']

assert(hm1.IsCompatible(hm2))
assert(np.allclose(hm1, hm2))

print('Checking healpy.read_map')
hm3 = hp.read_map('test_map.fits')
hm3 = coordinateutils.HealpixSkyMap(hm3)
assert(hm1.IsCompatible(hm3))
assert(np.allclose(hm1, hm3))

print('Checking healpy.write_map')
os.remove('test_map.fits')
hp.write_map('test_map.fits', np.asarray(hm3))
hm4 = coordinateutils.maputils.load_skymap_fits('test_map.fits')['T']
assert(hm1.IsCompatible(hm4))
assert(np.allclose(hm1, hm4))

os.remove('test_map.fits')
