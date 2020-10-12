#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
from spt3g import core, maps
from astropy.wcs import WCS
from astropy.io import fits

for dim in [300, 301]:

    arr = np.arange(dim * dim).reshape(dim, dim)
    res = core.G3Units.arcmin
    deg = core.G3Units.deg
    a0 = 20 * deg
    d0 = -50 * deg
    d1 = -51 * deg
    x0 = 120
    y0 = 180

    verbose = False
    error = ''

    fm0 = maps.FlatSkyMap(dim, dim, res, proj=maps.MapProjection.ProjZEA,
                          alpha_center=a0, delta_center=d0,
                          x_center=x0, y_center=y0)
    fm0b = maps.FlatSkyMap(dim, dim, res, proj=maps.MapProjection.ProjZEA,
                           alpha_center=a0, delta_center=d1,
                           x_center=x0, y_center=y0)

    print('Checking projections for dim {}'.format(dim))
    if verbose:
        print('-' * 80)
    for p in [0, 1, 2, 3, 4, 5, 6, 7, 9]:
        if verbose:
            print('Checking Proj{}'.format(p))

        proj = getattr(maps.MapProjection, 'Proj{}'.format(p))
        fm1 = maps.FlatSkyMap(arr, res, proj=proj,
                              alpha_center=a0, delta_center=d0,
                              x_center=x0, y_center=y0)
        assert(np.allclose(fm1.angle_to_xy(a0, d0), (x0, y0)))
        assert(np.allclose(fm1.xy_to_angle(x0, y0), (a0, d0)))
        assert(fm1.alpha_center == a0)
        assert(fm1.delta_center == d0)
        assert(fm1.x_center == x0)
        assert(fm1.y_center == y0)

        fmc = fm0.clone(False)
        fmc.proj = proj
        assert(fmc.compatible(fm1))

        fmcb = fm0b.clone(False)
        fmcb.proj = proj
        fmc.delta_center = d1
        assert(fmcb.compatible(fmc))
        assert(np.all(np.array(fmcb.angle_to_xy(a0, d0)) == np.array(fmc.angle_to_xy(a0, d0))))

        maps.fitsio.save_skymap_fits('test_map.fits', fm1, overwrite=True)
        fm2 = maps.fitsio.load_skymap_fits('test_map.fits')['T']

        try:
            assert(fm1.compatible(fm2))
        except AssertionError:
            for attr in ['alpha_center', 'delta_center', 'x_center', 'y_center']:
                print(attr, getattr(fm1, attr), getattr(fm2, attr))
            error += '\nProj{}: fits error'.format(p)
        assert(np.allclose(fm1, fm2))
        assert(np.allclose(fm1.wcs.wcs.crpix, fm2.wcs.wcs.crpix))
        assert(np.allclose(fm1.wcs.wcs.crval, fm2.wcs.wcs.crval))
        assert(np.allclose(fm1.wcs.wcs.cdelt, fm2.wcs.wcs.cdelt))
        assert(np.allclose(fm1.wcs.wcs.pc, fm2.wcs.wcs.pc))

        w = fm2.wcs
        hdr = maps.fitsio.create_wcs_header(fm2)
        if verbose:
            print(repr(hdr))
        pixs = np.array([[0, 0], [0, dim - 1], [dim // 2, dim // 2], [dim - 1, 0], [dim - 1, dim - 1]]).astype(float)
        ra, dec = maps.get_ra_dec_map(fm2)

        angs = []
        bad = False
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
                if verbose:
                    raise AssertionError
            except AssertionError:
                if not bad and not verbose:
                    print('ERROR: Proj{}'.format(p))
                    print('-' * 80)
                    print(repr(hdr))
                    print('-' * 80)
                bad = True
                print('pix', pix, 'pix2ang', [ra[idx] / deg, dec[idx] / deg],
                      'g3', g3_ang, 'wcs', wcs_ang, 'diff', np.abs(g3_ang - wcs_ang))
                error += '\nProj{}: xy_to_angle error'.format(p)
            angs.append(g3_ang * deg)

        for ang, pix in zip(angs, pixs):
            wcs_pix = np.asarray(w.all_world2pix(ang[0] / deg, ang[1] / deg, 0))
            g3_pix = np.asarray(fm2.angle_to_xy(*ang))
            try:
                assert(np.allclose(g3_pix, pix)) # round trip
                assert(np.allclose(wcs_pix, g3_pix))
                if verbose:
                    raise AssertionError
            except AssertionError:
                if not bad and not verbose:
                    print('ERROR: Proj{}'.format(p))
                    print('-' * 80)
                    print(repr(hdr))
                    print('-' * 80)
                bad = True
                print('ang', ang / deg, 'pix', pix, 'g3', g3_pix, 'wcs', wcs_pix,
                      'diff', np.abs(g3_pix - wcs_pix))
                error += '\nProj{}: angle_to_xy error'.format(p)

        if bad or verbose:
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
    assert(hm1.pol_conv == hm2['U'].pol_conv)
    assert(hm1.compatible(hm2['T']))
    assert(np.allclose(hm1, hm2['T']))
    assert(np.allclose(hm1, hm2['U']))

    print('Checking healpy.read_map')
    hm3 = hp.read_map('test_map.fits')
    hm3 = maps.HealpixSkyMap(hm3)
    assert(hm1.compatible(hm3))
    assert(np.allclose(hm1, hm3))

    print('Checking healpy.write_map')
    os.remove('test_map.fits')
    hp.write_map('test_map.fits', np.asarray(hm3))
    hm4 = maps.fitsio.load_skymap_fits('test_map.fits')['T']
    assert(hm1.compatible(hm4))
    assert(np.allclose(hm1, hm4))

finally:
    if os.path.exists('test_map.fits'):
        os.remove('test_map.fits')
