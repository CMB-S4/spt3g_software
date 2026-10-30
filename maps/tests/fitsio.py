#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
import gzip
import shutil
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

    # empty sparse maps
    sm = maps.HealpixSkyMap(64, pol_conv=maps.MapPolConv.IAU)
    maps.fitsio.save_skymap_fits('test_map.fits', sm, sm, sm, overwrite=True)
    out = maps.fitsio.load_skymap_fits('test_map.fits')
    assert(all(out[k].npix_nonzero == 0 and out[k].nside == 64 for k in 'TQU'))

    print('Checking chunked table IO')
    nrows = 1000
    tcols = [
        fits.Column(name='PIXEL', format='J', array=np.arange(nrows) * 3),
        fits.Column(name='A', format='D', array=rng.normal(size=nrows)),
        fits.Column(name='B', format='4E', array=rng.normal(size=(nrows, 4))),
    ]
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(tcols)]).writeto(
        'test_map.fits', overwrite=True)
    with fits.open('test_map.fits') as hl:
        rowdtype = hl[1].columns.dtype.newbyteorder('>')
        offset = hl.fileinfo(1)['datLoc']
        tdata = {c.name: np.array(hl[1].data[c.name]) for c in tcols}
    # chunks of 7 rows, so that the last chunk is partial
    chunk = 7 * rowdtype.itemsize
    for name, dtype in [('PIXEL', np.int64), ('A', float), ('B', float)]:
        col = maps.fitsio.read_table_column(
            'test_map.fits', offset, rowdtype, nrows, name, dtype, chunk=chunk
        )
        assert(col.dtype == dtype and np.array_equal(col, tdata[name].ravel()))
    with open('test_map.fits', 'r+b') as f:
        f.seek(offset)
        f.write(bytes(nrows * rowdtype.itemsize))
        for name in tdata:
            maps.fitsio.write_table_column(
                f, offset, rowdtype, name, tdata[name] * 2, ufac=2, chunk=chunk
            )
    with fits.open('test_map.fits') as hl:
        assert(all(np.allclose(hl[1].data[k], v) for k, v in tdata.items()))

    # compressed files are read through astropy
    maps.fitsio.save_skymap_fits('test_map.fits', sT, sU, sU, overwrite=True)
    with open('test_map.fits', 'rb') as fi, gzip.open('test_map.fits.gz', 'wb') as fo:
        shutil.copyfileobj(fi, fo)
    out = maps.fitsio.load_skymap_fits('test_map.fits.gz')
    for k, sm in zip('TQU', (sT, sU, sU)):
        assert(same_sparse(out[k], sm))

    default_chunk = maps.fitsio.chunk_size
    # chunks of 7 image rows, so that the last chunk is partial
    maps.fitsio.chunk_size = 7 * 8 * 500
    try:
        print('Checking chunked flat map writes')
        fmap = maps.FlatSkyMap(500, 300, U.arcmin, proj=maps.MapProjection.ProjZEA,
                               pol_conv=maps.MapPolConv.IAU, weighted=False)
        fmap[50:150, 100:400] = rng.normal(size=(100, 300))
        fdense = np.asarray(fmap.clone())
        fw = make_weights(fmap)
        for kwargs, atol in [
            (dict(W=fw), 1e-12),
            (dict(units='uK_CMB'), 1e-12),
            (dict(W=fw, compress=True), 0.1),
        ]:
            maps.fitsio.save_skymap_fits('test_map.fits', fmap, fmap, fmap, overwrite=True,
                                         **kwargs)
            # sparse input maps are not converted to dense maps
            assert(fmap.sparse and fw.TT.sparse)
            if kwargs.get('compress'):
                # all maps are compressed, following an empty primary HDU
                with fits.open('test_map.fits') as hl:
                    assert(hl[0].data is None and len(hl) == 10)
                    assert(all(isinstance(h, fits.CompImageHDU) for h in hl[1:]))
            out = maps.fitsio.load_skymap_fits('test_map.fits')
            # maps that were sparse when saved are loaded as sparse maps
            for k in 'TQU':
                assert(out[k].sparse)
                assert(np.allclose(np.asarray(out[k]), fdense, rtol=0, atol=atol))
            if 'W' in kwargs:
                assert(out['W'].TT.sparse)
                assert(np.allclose(np.asarray(out['W'].TT), 2 * fdense, rtol=0, atol=2 * atol))

        print('Checking chunked flat map reads')
        # irregular coverage, so that chunks have different nonzero spans
        yy, xx = np.mgrid[:300, :500]
        disc = (yy - 150) ** 2 + (xx - 250) ** 2 < 100 ** 2
        fdisc = maps.FlatSkyMap(np.where(disc, rng.normal(size=(300, 500)), 0.0),
                                U.arcmin, proj=maps.MapProjection.ProjZEA, weighted=False)
        fdisc.compact()
        ddisc = np.asarray(fdisc.clone())
        nzdisc = fdisc.npix_allocated
        for sparse in [True, False]:
            fdisc.sparse = sparse
            maps.fitsio.save_skymap_fits('test_map.fits', fdisc, overwrite=True)
            with fits.open('test_map.fits') as hl:
                # image stored with the axes swapped (DEC first), so that it is transposed
                thdr = hl[0].header.copy()
                thdr.update(WCS(thdr).swapaxes(0, 1).to_header())
                for k in ['X0', 'Y0', 'ALPHA0', 'DELTA0']:
                    thdr.remove(k, ignore_missing=True)
                fits.PrimaryHDU(np.ascontiguousarray(hl[0].data.T), header=thdr).writeto(
                    'test_map_t.fits', overwrite=True)
                # scaled image, read through astropy
                hl[0].scale('int32', bscale=1e-6, bzero=0.5)
                hl.writeto('test_map_s.fits', overwrite=True)
            with open('test_map.fits', 'rb') as fi, gzip.open('test_map.fits.gz', 'wb') as fo:
                shutil.copyfileobj(fi, fo)
            for fname, atol in [('test_map.fits', 0), ('test_map_t.fits', 0),
                                ('test_map.fits.gz', 0), ('test_map_s.fits', 1e-6)]:
                out = maps.fitsio.load_skymap_fits(fname)['T']
                # maps are loaded with the storage they had when saved
                assert(out.shape == fdisc.shape and out.sparse == sparse)
                assert(np.allclose(np.asarray(out.clone()), ddisc, rtol=0, atol=atol))
                if sparse and not atol:
                    # only the nonzero pixels are stored
                    assert(out.npix_allocated == nzdisc)

        print('Checking masks')
        # non-square flat mask, read and written in chunks of rows
        fmask = maps.G3SkyMapMask(fdisc, use_data=True)
        for kwargs in [dict(), dict(compress=True)]:
            maps.fitsio.save_skymapmask_fits('test_mask.fits', fmask, overwrite=True,
                                             **kwargs)
            with fits.open('test_mask.fits') as hl:
                # stored as unsigned bytes, with compressed images following an
                # empty primary HDU
                img = [h for h in hl if h.data is not None]
                assert(len(img) == 1 and img[0].data.dtype == np.uint8)
                assert(isinstance(img[0], fits.CompImageHDU) == bool(kwargs))
            out = maps.fitsio.load_skymapmask_fits('test_mask.fits')
            assert(out.compatible(fmask) and (out == fmask).all())
        # masks of transposed and scaled images, and masks stored as signed bytes
        # by earlier versions
        fits.PrimaryHDU(np.asarray(fmask), header=maps.fitsio.create_wcs_header(fdisc)
                        ).writeto('test_mask.fits', overwrite=True)
        for fname in ['test_map_t.fits', 'test_map_s.fits', 'test_mask.fits']:
            out = maps.fitsio.load_skymapmask_fits(fname)
            assert(out.compatible(fmask) and (out == fmask).all())
        # healpix mask, written in more than one chunk
        hmask = maps.G3SkyMapMask(sT, use_data=True)
        maps.fitsio.save_skymapmask_fits('test_mask.fits', hmask, overwrite=True)
        assert(fits.getheader('test_mask.fits', 1)['NAXIS2'] > maps.fitsio.chunk_size)
        out = maps.fitsio.load_skymapmask_fits('test_mask.fits')
        assert(out.compatible(hmask) and (out == hmask).all())
    finally:
        maps.fitsio.chunk_size = default_chunk

    hmask1 = maps.G3SkyMapMask(hm1, use_data=True)
    maps.fitsio.save_skymapmask_fits("test_mask.fits", hmask1, overwrite=True)
    hmask2 = maps.fitsio.load_skymapmask_fits("test_mask.fits")

    assert(hmask1.compatible(hmask2))
    assert((hmask1 == hmask2).all())

finally:
    for f in ['test_map.fits', 'test_map_t.fits', 'test_map_s.fits', 'test_map.fits.gz',
              'test_mask.fits']:
        if os.path.exists(f):
            os.remove(f)
