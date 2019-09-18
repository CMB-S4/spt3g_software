#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
from spt3g import core, coordinateutils
from astropy.wcs import WCS
from astropy.io import fits

arr = np.arange(300 * 300).reshape(300, 300)

for p in [0, 1, 2, 4, 5, 6, 7, 9]:
    print('Checking Proj{}'.format(p))

    proj = getattr(coordinateutils.MapProjection, 'Proj{}'.format(p))
    fm1 = coordinateutils.FlatSkyMap(arr, core.G3Units.arcmin, proj=proj,
                                     alpha_center=20 * core.G3Units.deg,
                                     delta_center=-50 * core.G3Units.deg)

    coordinateutils.maputils.save_skymap_fits('test_map.fits', fm1, overwrite=True)
    fm2 = coordinateutils.maputils.load_skymap_fits('test_map.fits')['T']

    assert(fm1.IsCompatible(fm2))
    assert(np.allclose(fm1, fm2))

    with fits.open('test_map.fits') as hdulist:
        pixs = np.asarray([[0, 0], [0, 300], [150, 150], [300, 0], [300, 300]], dtype=float)
        w = WCS(hdulist[-1].header)
        for pix in pixs:
            wcs_ang = np.asarray(w.wcs_pix2world([pix], 0))
            wcs_ang[wcs_ang > 180] -= 360
            g3_ang = np.asarray(fm1.xy_to_angle(*pix)) / core.G3Units.deg
            assert(np.allclose(wcs_ang, g3_ang))

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
