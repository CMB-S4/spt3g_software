#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
from spt3g import core, maps
from astropy.wcs import WCS
from astropy.io import fits


def assert_raises(exc, func, *args, **kwargs):
    """Check that the function call raises the given exception."""
    try:
        func(*args, **kwargs)
    except exc:
        return
    raise AssertionError("Expected {} from {}".format(exc.__name__, func.__name__))


def stored_units(fname):
    """Units label of the first map stored in the file."""
    with fits.open(fname) as hl:
        hdr = next(h.header for h in hl if h.data is not None)
    return hdr.get('BUNIT', hdr.get('TUNIT1'))


def make_weights(m, scale=2):
    """Polarized weights with each component equal to the map times the scale."""
    w = maps.G3SkyMapWeights(m)
    for k in maps.fitsio.weight_cols:
        setattr(w, k, m * scale)
    return w


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

        fmask1 = maps.G3SkyMapMask(fm1, use_data=True)
        maps.fitsio.save_skymapmask_fits("test_mask.fits", fmask1, overwrite=True)
        fmask2 = maps.fitsio.load_skymapmask_fits("test_mask.fits")

        assert(fmask1.compatible(fmask2))
        assert((fmask1 == fmask2).all())

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
    os.remove('test_mask.fits')
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
    hm3 = hp.read_map('test_map.fits', dtype=float)
    hm3 = maps.HealpixSkyMap(hm3)
    assert(hm1.compatible(hm3))
    assert(np.allclose(hm1, hm3))

    print('Checking units of healpy files')
    U = core.G3Units
    TU = core.G3TimestreamUnits
    # (units stored in the file, units argument, expected unit value or None for error)
    for column_units, units, factor in [
        (None, None, None),
        (None, 'mK_CMB', U.mK),
        (None, 'uK_CMB', U.uK),
        ('K_CMB', None, U.K),
        ('K_CMB', 'K_CMB', U.K),
        ('K_CMB', 'uK_CMB', None),
        ('Jy', 'Jy', U.Jy),
        ('Jy', 'mJy', None),
    ]:
        hp.write_map('test_map.fits', np.asarray(hm3), dtype=float, coord="C",
                     column_units=column_units, overwrite=True)
        load = maps.fitsio.load_skymap_fits
        if factor is None:
            assert_raises(ValueError, load, 'test_map.fits', units=units)
            continue
        hm4 = load('test_map.fits', units=units)['T']
        assert(hm1.compatible(hm4) and not hm4.weighted)
        assert(hm4.units == maps.fitsio.parse_fits_units(units or column_units)[0])
        assert(np.allclose(np.asarray(hm1) * factor, hm4, rtol=1e-12, atol=0))

    print('Checking units parsing')
    for ustr, utype, factor in [
        ('K_CMB', TU.Tcmb, U.K),
        ('KCMB', TU.Tcmb, U.K),
        ('uK_CMB', TU.Tcmb, U.uK),
        ('muK', TU.Tcmb, U.uK),
        ('mK_RJ', TU.Trj, U.mK),
        ('K_CMB^2', TU.Tcmb, U.K ** 2),
        ('(uK_CMB)^2', TU.Tcmb, U.uK ** 2),
        ('mJy', TU.FluxDensity, U.mJy),
        ('MJy/sr', TU.FluxDensity, U.MJy / U.sr),
        ('Tcmb', TU.Tcmb, 1.0),
        ('None', getattr(TU, 'None'), 1.0),
    ]:
        out = maps.fitsio.parse_fits_units(ustr)
        assert(out[0] == utype and np.isclose(out[1], factor)), (ustr, out)
    for ustr in ['furlongs', 'K/sr', 'Jy_CMB', 'K_CMB^x']:
        assert_raises(ValueError, maps.fitsio.parse_fits_units, ustr)

    print('Checking weights and units round trip')
    for m1 in [
        maps.HealpixSkyMap(np.arange(12 * 16 * 16) + 1.0, pol_conv=maps.MapPolConv.IAU,
                           weighted=False),
        maps.FlatSkyMap(np.arange(300.).reshape(15, 20) + 1, U.arcmin,
                        proj=maps.MapProjection.ProjZEA,
                        pol_conv=maps.MapPolConv.IAU, weighted=False),
    ]:
        m1w = m1.clone()
        m1w.weighted = True
        w = make_weights(m1)

        # maps saved without weights are stored as unweighted maps
        for m, wt, label, weighted in [
            (m1, None, 'mK_CMB', False),
            (m1, w, 'mK_CMB', False),
            (m1w, None, 'mK_CMB', False),
            (m1w, w, 'Tcmb', True),
        ]:
            maps.fitsio.save_skymap_fits('test_map.fits', m, m, m, W=wt, overwrite=True)
            assert(stored_units('test_map.fits') == label)
            out = maps.fitsio.load_skymap_fits('test_map.fits')
            assert(('W' in out) == (wt is not None))
            for k in 'TQU':
                assert(out[k].units == TU.Tcmb and out[k].weighted == weighted)
                assert(np.allclose(out[k], m1))
            if wt is not None:
                assert(out['W'].polarized)
                for k in maps.fitsio.weight_cols:
                    assert(np.allclose(getattr(out['W'], k), getattr(w, k)))
        # the weighted flag does not depend on whether the weights are loaded
        out = maps.fitsio.load_skymap_fits('test_map.fits', keys=['T', 'Q', 'U'])
        assert('W' not in out and all(out[k].weighted for k in 'TQU'))

        print('Checking output units for {}'.format(type(m1).__name__))
        m2 = m1.clone()
        m2.units = TU.FluxDensity
        m2 *= U.mJy
        # (map, units argument, stored label, unit value)
        for m, ustr, label, uval in [
            (m2, None, 'Jy', U.Jy),
            (m1, 'uK_CMB', 'uK_CMB', U.uK),
            (m1w, 'uK_CMB', 'uK_CMB', U.uK),
            (m2, 'mJy', 'mJy', U.mJy),
        ]:
            maps.fitsio.save_skymap_fits('test_map.fits', m, overwrite=True, units=ustr)
            assert(stored_units('test_map.fits') == label)
            with fits.open('test_map.fits') as hl:
                raw = next(h.data for h in hl if h.data is not None)
            if isinstance(m, maps.HealpixSkyMap):
                raw = raw['TEMPERATURE']
            assert(np.allclose(raw, np.asarray(m) / uval, rtol=1e-12, atol=0))
            out = maps.fitsio.load_skymap_fits('test_map.fits')['T']
            assert(out.units == m.units and not out.weighted)
            assert(np.allclose(out, m, rtol=1e-12, atol=0))
        for m, kwargs in [
            (m1, dict(units='mJy')),
            (m1, dict(units='uK_RJ')),
            (m1w, dict(units='uK_CMB', W=w)),
        ]:
            assert_raises(ValueError, maps.fitsio.save_skymap_fits, 'test_map.fits',
                          m, m, m, overwrite=True, **kwargs)

    print('Checking sparse healpix maps')

    def same_sparse(a, b):
        pa, va = (np.asarray(x) for x in a.nonzero_pixels())
        pb, vb = (np.asarray(x) for x in b.nonzero_pixels())
        ia, ib = np.argsort(pa), np.argsort(pb)
        return np.array_equal(pa[ia], pb[ib]) and np.allclose(va[ia], vb[ib])

    rng = np.random.default_rng(42)
    spix = rng.choice(12 * 64 * 64, 2000, replace=False)
    extra = int(np.setdiff1d(np.arange(12 * 64 * 64), spix)[0])
    for nested in [False, True]:
        smaps = []
        for scale in [1, 2, 3]:
            sm = maps.HealpixSkyMap(64, nested=nested, pol_conv=maps.MapPolConv.IAU)
            for p in spix:
                sm[int(p)] = scale * (1.0 + p)
            smaps.append(sm)
        sT, sQ, sU = smaps
        # store U with a different sparse layout, which lists the same pixels in a
        # different order, to exercise the per-map reordering
        if sU.ringsparse:
            sU.indexedsparse = True
        else:
            sU.ringsparse = True
        assert(not np.array_equal(*(np.asarray(m.nonzero_pixels()[0]) for m in (sT, sU))))
        sw = make_weights(sT)

        maps.fitsio.save_skymap_fits('test_map.fits', sT, sQ, sU, W=sw, overwrite=True)
        filepix = fits.getdata('test_map.fits', 1)['PIXEL']
        assert(len(filepix) == len(spix) and np.all(np.diff(filepix) > 0))
        out = maps.fitsio.load_skymap_fits('test_map.fits')
        for k, sm in zip('TQU', smaps):
            assert(not out[k].dense and out[k].nested == nested)
            assert(same_sparse(out[k], sm))
        for k in maps.fitsio.weight_cols:
            assert(same_sparse(getattr(out['W'], k), getattr(sw, k)))

        # maps with different nonzero pixels than T
        sQ[extra] = 1.0
        assert_raises(ValueError, maps.fitsio.save_skymap_fits, 'test_map.fits',
                      sT, sQ, sU, overwrite=True)

    hmask1 = maps.G3SkyMapMask(hm1, use_data=True)
    maps.fitsio.save_skymapmask_fits("test_mask.fits", hmask1, overwrite=True)
    hmask2 = maps.fitsio.load_skymapmask_fits("test_mask.fits")

    assert(hmask1.compatible(hmask2))
    assert((hmask1 == hmask2).all())

finally:
    if os.path.exists('test_map.fits'):
        os.remove('test_map.fits')
    if os.path.exists('test_mask.fits'):
        os.remove('test_mask.fits')
