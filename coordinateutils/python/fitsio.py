from spt3g import core
from spt3g.coordinateutils import HealpixSkyMap, FlatSkyMap, G3SkyMapWeights
from spt3g.coordinateutils import MapPolType, WeightType, MapCoordReference, MapProjection

import numpy as np
import os

__all__ = [
    'load_skymap_fits',
    'save_skymap_fits',
    'SaveMapFrame',
    'load_proj_dict',
    'create_wcs_header',
    'parse_wcs_header',
]


@core.usefulfunc
def load_skymap_fits(filename, hdu=None):
    """
    Load a fits file containing a sky map.

    Arguments
    ---------
    filename : str
        Path to fits file
    hdu : int
        If supplied, the data are extract from the given HDU index.

    Returns
    -------
    a dictionary of maps keyed with e.g. 'T', 'Q', 'U' and 'W'.
    """

    import astropy.io.fits
    assert(os.path.exists(filename))

    map_type = None
    weight_type = None
    map_opts = {}
    output = {}

    # defaults for missing header entries
    polcconv = 'COSMO'
    units = 'Tcmb'
    coord_ref = 'Equatorial'
    proj = 'Proj5'
    is_weighted = False
    alpha_center = None
    delta_center = None
    res = None
    xres = None

    with astropy.io.fits.open(filename) as hdulist:
        for hidx, H in enumerate(hdulist):

            hdr = H.header
            mtype = hdr.get('MAPTYPE', None)

            if mtype == 'FLAT' or 'PROJ' in hdr or 'WCSAXES' in hdr or 'CTYPE1' in hdr:
                # flat map
                if not map_type:
                    map_type = 'flat'
                elif map_type != 'flat':
                    raise ValueError(
                        "Expected a {} sky map in HDU {}, found a flat map".format(
                            map_type, hidx
                        )
                    )
            elif mtype == 'HEALPIX' or hdr.get('PIXTYPE', None) == 'HEALPIX' or 'NSIDE' in hdr:
                # healpix map
                if not map_type:
                    map_type = 'healpix'
                elif map_type != 'healpix':
                    raise ValueError(
                        "Expected a {} sky map in HDU {}, found a healpix map".format(
                            map_type, hidx
                        )
                    )

            polcconv = hdr.get('POLCCONV', polcconv)
            if 'TUNIT' in hdr:
                udict = {'k_cmb': 'Tcmb'}
                units = udict.get(hdr['TUNIT'].lower(), units)
            else:
                units = hdr.get('UNITS', units)
            if 'COORDSYS' in hdr:
                cdict = {'C': 'Equatorial', 'G': 'Galactic', 'L': 'Local'}
                coord_ref = cdict.get(hdr['COORDSYS'], coord_ref)
            else:
                coord_ref = hdr.get('COORDREF', coord_ref)
            overflow = hdr.get('OVERFLOW', 0)

            map_opts.update(
                units=getattr(core.G3TimestreamUnits, units),
                coord_ref=getattr(MapCoordReference, coord_ref),
            )

            if map_type == 'flat':
                map_opts.update(parse_wcs_header(hdr))

            elif map_type == 'healpix':
                nside = hdr['NSIDE']
                map_opts.update(is_nested=(hdr['ORDERING'] == 'nest'))

            # primary HDU
            if H.data is None:
                continue

            # extracting a particular HDU
            if map_type is not None and hdu is not None and hidx != hdu:
                continue

            if map_type == 'flat':
                data = np.array(H.data, dtype=float)
                if map_opts.pop('transpose', False):
                    data = np.array(data.T)

            # map type must be known if the HDU contains map data
            if map_type not in ['flat', 'healpix']:
                raise ValueError("Unknown map type in HDU {}".format(hidx))

            if map_type == 'flat' and hdr.get('ISWEIGHT', None):
                # flat map weights
                assert('T' in output)
                if weight_type is None:
                    weight_type = 'Wpol' if ('Q' in output and 'U' in output) else 'Wunpol'
                    weight_type = getattr(WeightType, hdr.get('WGTTYPE', weight_type))
                weight_map = output.setdefault(
                    'W', G3SkyMapWeights(output['T'], weight_type)
                )
                fm = FlatSkyMap(data, **map_opts)
                fm.overflow = overflow
                setattr(weight_map, hdr['WTYPE'], fm)
                del data

            elif map_type == 'flat' and not hdr.get('ISWEIGHT', None):
                # flat map data
                ptype = hdr.get('POLTYPE', 'T')
                pol_type = getattr(MapPolType, ptype)
                if pol_type == MapPolType.U and polcconv == 'COSMO':
                    data *= -1
                fm = FlatSkyMap(data, pol_type=pol_type, **map_opts)
                fm.overflow = overflow
                output[ptype] = fm
                del data

            elif map_type == 'healpix':
                # healpix map data

                col_dict = {
                    'TEMPERATURE': 'T',
                    'Q_POLARISATION': 'Q',
                    'U_POLARISATION': 'U',
                    'I_STOKES': 'T',
                    'Q_STOKES': 'Q',
                    'U_STOKES': 'U',
                    'I': 'T',
                    'II': 'TT',
                    'IQ': 'TQ',
                    'IU': 'TU',
                }

                unit_dict = {
                    'k_cmb': 'Tcmb',
                }

                pix = None

                partial = hdr.get('INDXSCHM') == 'EXPLICIT' or hdr.get('OBJECT') == 'PARTIAL'

                for cidx, hcol in enumerate(H.data.names):
                    col = col_dict.get(hcol, hcol)
                    data = np.array(H.data[hcol], dtype=float).ravel()

                    if col == 'PIXEL' or (partial and cidx == 0):
                        pix = np.array(data, dtype=int).ravel()
                        del data
                        continue

                    units = unit_dict.get(hdr.get('TUNIT{:d}'.format(cidx + 1), units), units)
                    overflow = hdr.get('TOFLW{:d}'.format(cidx + 1), overflow)

                    is_weight = hdr.get(
                        'TISWGT{:d}'.format(cidx + 1),
                        col in ['TT', 'TQ', 'TU', 'QQ', 'QU', 'UU'],
                    )
                    if is_weight:
                        assert('T' in output)
                        if weight_type is None:
                            weight_type = 'Wpol' if ('Q' in output and 'U' in output) else 'Wunpol'
                            weight_type = getattr(WeightType, hdr.get('WGTTYPE', weight_type))
                        weight_map = output.setdefault(
                            'W', G3SkyMapWeights(output['T'], weight_type)
                        )

                        mdata = (pix, data, nside) if pix is not None else data
                        hm = HealpixSkyMap(mdata, **map_opts)
                        hm.overflow = overflow

                        setattr(weight_map, col, hm)

                    else:
                        pol_type = getattr(MapPolType, col)
                        if pol_type == MapPolType.U and polcconv == 'COSMO':
                            data *= -1
                        mdata = (pix, data, nside) if pix is not None else data

                        hm = HealpixSkyMap(mdata, pol_type=pol_type, **map_opts)
                        output[col] = hm

                    del mdata, data

                del pix
                break

            del H.data

    for k, m in output.items():
        if k == 'W':
            continue
        m.is_weighted = 'W' in output

    return output


@core.usefulfunc
def load_proj_dict(inverse=False):
    """
    Dictionary that associates projection names to their three-letter codes in WCS.
    """

    projdict = {
        'Proj0': 'SFL',
        'Proj1': 'CAR',
        'Proj2': 'SIN',
        'Proj3': '!!!',
        'Proj4': 'STG',
        'Proj5': 'ZEA',
        'Proj6': 'TAN',
        'Proj7': 'CEA',
        'Proj8': '!!!',
        'Proj9': 'CAR',
        'ProjSansonFlamsteed': 'SFL',
        'ProjSFL': 'SFL',
        'ProjPlateCarree': 'CAR',
        'ProjCAR': 'CAR',
        'ProjOrthographic': 'SIN',
        'ProjSIN': 'SIN',
        'ProjStereographic': 'STG',
        'ProjSTG': 'STG',
        'ProjLambertAzimuthalEqualArea': 'ZEA',
        'ProjZEA': 'ZEA',
        'ProjGnomonic': 'TAN',
        'ProjTAN': 'TAN',
        'ProjCylindricalEqualArea': 'CEA',
        'ProjCEA': 'CEA',
        'ProjBICEP': 'CAR',
        'ProjHealpix': 'HPX',
    }

    if inverse:
        projdict_inverse = {}
        for k, v in projdict.items():
            projdict_inverse.setdefault(v, []).append(k)
        return projdict_inverse

    return projdict


def get_wcs(skymap):
    """
    Creates WCS (world coordinate system) information from information in the input FlatSkyMap object.
    """

    if hasattr(skymap, '_wcs'):
        return skymap._wcs

    projdict = load_proj_dict()
    proj_abbr = projdict[str(skymap.proj)]

    from astropy.wcs import WCS

    w = WCS(naxis=2)

    if skymap.coord_ref == MapCoordReference.Equatorial:
        w.wcs.ctype = ['RA---{}'.format(proj_abbr), 'DEC--{}'.format(proj_abbr)]
        w.wcs.radesys = 'FK5'
        w.wcs.equinox = 2000
    elif skymap.coord_ref == MapCoordReference.Galactic:
        w.wcs.ctype = ['GLON-{}'.format(proj_abbr), 'GLAT-{}'.format(proj_abbr)]
    elif skymap.coord_ref == MapCoordReference.Local:
        w.wcs.ctype = ['RA---{}'.format(proj_abbr), 'DEC--{}'.format(proj_abbr)]
        w.wcs.radesys = 'ALTAZ'

    w.wcs.cdelt = [
        -skymap.x_res / core.G3Units.deg,
        skymap.res / core.G3Units.deg,
    ]
    w.wcs.cunit = ['deg', 'deg']

    crpix = [
        skymap.shape[1] - skymap.x_center + 1,
        skymap.shape[0] - skymap.y_center + 1
    ]
    crval = [
        skymap.alpha_center / core.G3Units.deg,
        skymap.delta_center / core.G3Units.deg,
    ]
    if proj_abbr in ['CEA']:
        v = skymap.delta_center / core.G3Units.rad
        w.wcs.set_pv([(2, 1, np.cos(v) ** 2)])
        crval[1] = 0.0
    w.wcs.crpix = crpix
    w.wcs.crval = crval

    skymap._wcs = w
    return w

# set property
setattr(FlatSkyMap, 'wcs', property(get_wcs))


@core.usefulfunc
def create_wcs_header(skymap):
    return skymap.wcs.to_header()


@core.usefulfunc
def parse_wcs_header(header):
    """
    Extract flat sky map keywords from a WCS fits header.  Raises an error if
    the header is malformed.
    """

    from astropy.wcs import WCS
    w = WCS(header)

    # parse projection
    ctype = w.wcs.ctype
    wcsproj = ctype[0][-3:]
    proj = header.get('PROJ', None)
    projdict = load_proj_dict(inverse=True)
    if wcsproj not in projdict or wcsproj == '!!!':
        raise ValueError('Unknown WCS projection {}'.format(wcsproj))
    projopts = projdict[wcsproj]
    if proj is None:
        proj = projopts[0]
    elif proj not in projopts:
        raise ValueError('PROJ keyword inconsistent with WCS coordinate type')

    # parse coordinate system
    coord_ref = None
    if ctype[0].startswith('RA-') or ctype[0].startswith('DEC-'):
        if w.wcs.radesys == 'FK5':
            coord_ref = MapCoordReference.Equatorial
        elif w.wcs.radesys == 'ALTAZ':
            coord_ref = MapCoordReference.Local
    elif ctype[0].startswith('GLON-') or ctype[0].startswith('GLAT-'):
        coord_ref = MapCoordReference.Galactic

    # parse coordinate order
    if ctype[0].startswith('RA-') or ctype[0].startswith('GLON-'):
        order = [0, 1]
        transpose = False
    elif ctype[0].startswith('DEC-') or ctype[0].startswith('GLAT-'):
        order = [1, 0]
        transpose = True
    else:
        raise ValueError('Unable to determine WCS coordinate axes')

    # parse resolution
    cdelt = w.wcs.get_cdelt() * np.diag(w.wcs.get_pc())
    x_res = np.abs(cdelt[order[0]]) * core.G3Units.deg
    res = np.abs(cdelt[order[1]]) * core.G3Units.deg

    # parse map center
    crval = w.wcs.crval
    alpha_center = crval[order[0]] * core.G3Units.deg
    delta_center = crval[order[1]] * core.G3Units.deg

    crpix = w.wcs.crpix
    x_center = header['NAXIS{}'.format(order[0] + 1)] - crpix[order[0]] + 1
    y_center = header['NAXIS{}'.format(order[1] + 1)] - crpix[order[1]] + 1

    if wcsproj in ['CEA']:
        for i, m, v in w.wcs.get_pv():
            if i == order[1] + 1 and m == 1:
                delta_center = np.arccos(np.sqrt(v)) * core.G3Units.rad

    # construct arguments
    map_opts = dict(
        proj=getattr(MapProjection, proj),
        res=res,
        x_res=x_res,
        alpha_center=alpha_center,
        delta_center=delta_center,
        x_center=x_center,
        y_center=y_center,
        transpose=transpose,
    )
    if coord_ref is not None:
        map_opts.update(coord_ref=coord_ref)

    return map_opts


@core.usefulfunc
def save_skymap_fits(filename, T, Q=None, U=None, W=None, overwrite=False,
                     compress=True):
    """
    Save G3 map objects to a fits file.

    `FlatSkyMap` objects are stored in a series of (optionally compressed)
    `ImageHDU` entries, in which each HDU contains the projection information in
    its header in standard WCS format, along with the image data for a single
    map (one of the Stokes maps or a weight map component).

    `HealpixSkyMap` objects are stored in a `BinTableHDU` extension, which
    contains the necessary header information for compatiblity with healpix map
    readers (e.g. `healpix.read_map`), and a single table with one column per
    Stokes map or weight map component.  Sparse maps are stored as cut-sky
    pixel-indexed tables, while dense maps are stored with implicit indexing
    over all pixels.  The former produces output that is equivalent to using
    `healpy.write_map` with the `partial=True` option.

    Arguments
    ---------
    filename : str
        Path to output file.  Must not exist, unless overwrite is True.
    T[, Q, U] : FlatSkyMap or HealpixSkyMap
        Maps to save
    W : G3SkyMapWeights
        Weights to save with the maps
    overwrite : bool
        If True, any existing file with the same name will be ovewritten.
    compress : bool
        If True, and if input maps are FlatSkyMap objects, store these in a
        series of compressed image HDUs, one per map.  Otherwise, store input
        maps in a series of standard ImageHDUs, which are readable with older
        FITS readers (e.g. idlastro).
    """

    import astropy.io.fits

    if isinstance(T, (FlatSkyMap, HealpixSkyMap)):
        flat = isinstance(T, FlatSkyMap)
    else:
        raise TypeError("Input map must be a FlatSkyMap or HealpixSkyMap instance")

    if Q is not None:
        assert(U is not None)
        assert(T.IsCompatible(Q))
        assert(T.IsCompatible(U))
        pol = True
        maps = [T, Q, U]
        names = ['T', 'Q', 'U']
    else:
        assert(U is None)
        pol = False
        maps = [T]
        names = ['T']

    if flat:
        header = create_wcs_header(T)

        header['PROJ'] = str(T.proj)

        bitpix = {
            np.dtype(np.int16): 16,
            np.dtype(np.int32): 32,
            np.dtype(np.float32): -32,
            np.dtype(np.float64): -64,
        }
        header['BITPIX'] = bitpix[np.asarray(T).dtype]
        header['NAXIS'] = 2

    else:
        header = astropy.io.fits.Header()

        header['PIXTYPE'] = 'HEALPIX'
        header['ORDERING'] = 'NEST' if T.nested else 'RING'
        header['NSIDE'] = T.nside
        if T.dense:
            header['INDXSCHM'] = 'IMPLICIT'
            header['OBJECT'] = 'FULLSKY'
        else:
            header['INDXSCHM'] = 'EXPLICIT'
            header['OBJECT'] = 'PARTIAL'

        cdict = {'Equatorial': 'C', 'Galactic': 'G', 'Local': 'L'}
        header['COORDSYS'] = cdict[str(T.coord_ref)]
        if header['COORDSYS'] == 'C':
            header['RADESYS'] = 'FK5'
            header['EQUINOX'] = 2000.0

        conv = {
            np.dtype(np.int16): 'I',
            np.dtype(np.int32): 'J',
            np.dtype(np.int64): 'K',
            np.dtype(np.float32): 'E',
            np.dtype(np.float64): 'D',
        }
        pix = None

    header['MAPTYPE'] = 'FLAT' if flat else 'HEALPIX'
    header['COORDREF'] = str(T.coord_ref)
    header['POLCCONV'] = 'IAU'
    header['UNITS'] = str(T.units)
    header['WEIGHTED'] = T.is_weighted
    if W is not None:
        header['WGTTYPE'] = str(W.weight_type)

    hdulist = astropy.io.fits.HDUList()
    cols = astropy.io.fits.ColDefs([])

    try:
        for m, name in zip(maps, names):
            if flat:
                if compress:
                    hdu = astropy.io.fits.CompImageHDU(np.asarray(m), header=header)
                else:
                    hdu = astropy.io.fits.ImageHDU(np.asarray(m), header=header)
                hdu.header['ISWEIGHT'] = False
                hdu.header['POLTYPE'] = name
                hdu.header['OVERFLOW'] = m.overflow
                hdulist.append(hdu)
                del hdu
            else:
                if not T.dense:
                    pix1, data = m.nonzero_pixels()
                    pix1 = np.asarray(pix1)
                    data  = np.asarray(data)

                    if pix is None:
                        pix = pix1

                        fmt = conv.get(np.min_scalar_type(-pix.max()), 'K')
                        col = astropy.io.fits.Column(
                            name='PIXEL', format=fmt, array=pix, unit=None
                        )
                        cols.add_col(col)
                        del col
                    else:
                        # XXX different sparse maps can have different nonzero pixels
                        assert(len(pix) == len(pix1) and (pix == pix1).all())
                    del pix1

                else:
                    data = np.asarray(m)

                fmt = conv.get(data.dtype, 'D')
                col = astropy.io.fits.Column(
                    name=name, format=fmt, array=data, unit=str(m.units)
                )
                cols.add_col(col)
                del col
                del data

                idx = len(cols)
                header['TISWGT{:d}'.format(idx)] = False
                header['TOFLW{:d}'.format(idx)] = m.overflow

        if W is not None:
            if pol:
                assert(W.weight_type == WeightType.Wpol)
                wnames = ['TT', 'TQ', 'TU', 'QQ', 'QU', 'UU']
            else:
                assert(W.weight_type == WeightType.Wunpol)
                wnames = ['TT']

            for wt in wnames:
                m = getattr(W, wt)
                assert T.IsCompatible(m), 'Weight {} is not compatible with T'.format(wt)

                if flat:
                    if compress:
                        hdu = astropy.io.fits.CompImageHDU(np.asarray(m), header=header)
                    else:
                        hdu = astropy.io.fits.ImageHDU(np.asarray(m), header=header)
                    hdu.header['ISWEIGHT'] = True
                    hdu.header['WTYPE'] = wt
                    hdu.header['OVERFLOW'] = m.overflow
                    hdulist.append(hdu)
                    del hdu
                else:
                    if not T.dense:
                        pix1, data = m.nonzero_pixels()
                        pix1 = np.asarray(pix1)
                        data = np.asarray(data)
                        assert(len(pix) == len(pix1) and (pix == pix1).all())
                        del pix1
                    else:
                        data = np.asarray(m)

                    fmt = conv.get(data.dtype, 'D')
                    col = astropy.io.fits.Column(
                        name=wt, format=fmt, array=data, unit=str(m.units)
                    )
                    cols.add_col(col)
                    del col
                    del data

                    idx = len(cols)
                    header['TISWGT{:d}'.format(idx)] = False
                    header['TOFLW{:d}'.format(idx)] = m.overflow

        if not flat:
            del pix
            hdu = astropy.io.fits.BinTableHDU.from_columns(cols, header=header)
            hdulist.append(hdu)
            del hdu

        if overwrite:
            if os.path.exists(filename):
                os.remove(filename)
        else:
            assert(not os.path.exists(filename))

        hdulist.writeto(filename)

    finally:
        for col in cols.names:
            del cols[col].array
            cols.del_col(col)
        del cols
        for hdu in hdulist:
            del hdu.data
        del hdulist


@core.indexmod
def SaveMapFrame(frame, map_id, output_file, overwrite=False):
    """
    Save the map with Id map_id into output_file.
    """

    if frame.type != core.G3FrameType.Map:
        return

    if frame['Id'] != map_id:
        return

    T = frame['T']
    Q = frame.get('Q', None)
    U = frame.get('U', None)
    W = frame.get('Wpol', frame.get('Wunpol', None))

    save_skymap_fits(
        output_file,
        T=T,
        Q=Q,
        U=U,
        W=W,
        overwrite=overwrite,
    )
