#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
from spt3g import core, maps
from astropy.wcs import WCS
from astropy.io import fits

dim = 300
arr = np.arange(dim * dim).reshape(dim, dim)
res = core.G3Units.arcmin
deg = core.G3Units.deg
a0 = 20 * deg
d0 = -50 * deg
x0 = 120
y0 = 180

verbose = False
error = ''

print('Checking projections')
if verbose:
    print('-' * 80)
for p in [0, 1, 2, 4, 5, 6, 7, 9]:
    if verbose:
        print('Checking Proj{}'.format(p))

    proj = getattr(maps.MapProjection, 'Proj{}'.format(p))
    fm1 = maps.FlatSkyMap(arr, res, proj=proj,
                                     alpha_center=a0, delta_center=d0,
                                     x_center=x0, y_center=y0)

    maps.fitsio.save_skymap_fits('test_map.fits', fm1, overwrite=True)
    fm2 = maps.fitsio.load_skymap_fits('test_map.fits')['T']

    try:
        assert(fm1.IsCompatible(fm2))
    except AssertionError:
        for attr in ['alpha_center', 'delta_center', 'x_center', 'y_center']:
            print(attr, getattr(fm1, attr), getattr(fm2, attr))
        error += '\nProj{}: fits error'.format(p)
    assert(np.allclose(fm1, fm2))
    assert(fm1.wcs.to_header() == fm2.wcs.to_header())

    w = fm2.wcs
    if verbose:
        print(repr(w.to_header()))
    pixs = np.array([[0, 0], [0, 1], [0.5, 0.5], [1, 0], [1, 1]]) * dim
    ra, dec = maps.get_ra_dec_map(fm2)

    angs = []
    for pix in pixs:
        wcs_ang = np.asarray(w.all_pix2world(pix[0], pix[1], 0))
        wcs_ang[wcs_ang > 180] -= 360
        g3_ang = np.asarray(fm2.xy_to_angle(*pix)) / deg
        idx = (int(pix[1]), int(pix[0]))
        try:
            try:
                assert(np.allclose(g3_ang, [ra[idx] / deg, dec[idx] / deg]))
            except IndexError:
                pass
            assert(np.allclose(wcs_ang, g3_ang))
        except AssertionError:
            print(pix, g3_ang, wcs_ang, np.abs(g3_ang - wcs_ang))
            error += '\nProj{}: xy_to_angle error'.format(p)
        angs.append(g3_ang * deg)

    for ang, pix in zip(angs, pixs):
        wcs_pix = np.asarray(w.all_world2pix(ang[0] / deg, ang[1] / deg, 0))
        g3_pix = np.asarray(fm2.angle_to_xy(*ang))
        try:
            assert(np.allclose(g3_pix, pix)) # round trip
            assert(np.allclose(wcs_pix, g3_pix))
        except AssertionError:
            print(ang / deg, g3_pix, wcs_pix, np.abs(g3_pix - wcs_pix))
            error += '\nProj{}: angle_to_xy error'.format(p)

    if verbose:
        print('-' * 80)

if error:
    os.remove('test_map.fits')
    raise RuntimeError('Projection errors:{}'.format(error))

try:
    print('Checking Healpix')
    hm1 = maps.HealpixSkyMap(np.arange(12 * 64 * 64), pol_conv=maps.MapPolConv.COSMO)
    maps.fitsio.save_skymap_fits('test_map.fits', hm1, hm1, hm1, overwrite=True)
    hm2 = maps.fitsio.load_skymap_fits('test_map.fits')

    assert(hm2['T'].pol_conv == hm2['U'].pol_conv)
    assert(hm2['U'].pol_conv == maps.MapPolConv.IAU)
    assert(hm1.IsCompatible(hm2['T']))
    assert(np.allclose(hm1, hm2['T']))
    assert(np.allclose(hm1, -hm2['U']))

    print('Checking healpy.read_map')
    hm3 = hp.read_map('test_map.fits')
    hm3 = maps.HealpixSkyMap(hm3)
    assert(hm1.IsCompatible(hm3))
    assert(np.allclose(hm1, hm3))

    print('Checking healpy.write_map')
    os.remove('test_map.fits')
    hp.write_map('test_map.fits', np.asarray(hm3))
    hm4 = maps.fitsio.load_skymap_fits('test_map.fits')['T']
    assert(hm1.IsCompatible(hm4))
    assert(np.allclose(hm1, hm4))

finally:
    if os.path.exists('test_map.fits'):
        os.remove('test_map.fits')
