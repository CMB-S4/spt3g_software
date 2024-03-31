from spt3g import core
from spt3g.maps import HealpixSkyMap, FlatSkyMap, G3SkyMapWeights
from spt3g.maps import MapPolType, MapPolConv, MapCoordReference, MapProjection

import numpy as np
import os
import warnings

__all__ = [
    'load_skymap_fits',
    'save_skymap_fits',
    'SaveMapFrame',
]


@core.usefulfunc
def load_skymap_fits(filename, hdu=None, keys=None, memmap=False, apply_units=False):
    """
    Load a fits file containing a sky map.

    Arguments
    ---------
    filename : str
        Path to fits file
    hdu : int, optional
        If supplied, the data are extract from the given HDU index.
    keys : list of strings, optional
        If supplied, return only these keys in the output dictionary.
        Options are: T, Q, U, W.
    memmap : bool, optional
        Argument passed to astropy.io.fits.open. If True, the map is not read
        into memory, but only the required pixels are read when needed. Default:
        False.
    apply_units : bool, optional
        If True, and input maps have known units, multiply by the appropriate
        conversion factor to return maps in G3Units.

    Returns
    -------
    a dictionary of maps keyed with e.g. 'T', 'Q', 'U' and 'W'.
    """

    import astropy.io.fits
    assert(os.path.exists(filename))

    map_type = None
    pol = None
    map_opts = {}
    output = {}

    # defaults for missing header entries
    pol_conv = None
    units = 'Tcmb'
    coord_ref = 'Equatorial'
    proj = 'Proj5'
    weighted = False
    alpha_center = None
    delta_center = None
    res = None
    xres = None

    unit_dict = {
        'k_cmb': ('Tcmb', core.G3Units.K),
        'kcmb': ('Tcmb', core.G3Units.K),
        'kcmb^2': ('Tcmb', core.G3Units.K ** 2),
    }

    if keys is None:
        keys = ['T', 'Q', 'U', 'W']

    with astropy.io.fits.open(filename, memmap=memmap) as hdulist:
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
                # expect that flat sky maps on disk are probably in IAU
                if not pol_conv:
                    pol_conv = 'IAU'
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
                # expect that healpix maps on disk are probably in COSMO
                if not pol_conv:
                    pol_conv = 'COSMO'

            if 'POLAR' in hdr:
                pdict = {'T': True, 'F': False}
                pol = pdict.get(hdr['POLAR'], hdr['POLAR'])
            if pol and map_type is not None:
                if hdr.get('POLCCONV', '').upper() not in ['IAU', 'COSMO']:
                    core.log_warn('Polarization convention not set, assuming %s' % pol_conv,
                                  unit='load_skymap_fits')
                else:
                    pol_conv = hdr['POLCCONV'].upper()
            uconv = None
            if 'TUNIT' in hdr:
                u = hdr['TUNIT'].lower()
                if u in unit_dict:
                    units, uconv = unit_dict[u]
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
                pol_conv=getattr(MapPolConv, pol_conv) if pol_conv else None,
            )

            if map_type == 'flat':
                map_opts.update(flat_pol=hdr.get('FLATPOL', False))
                map_opts.update(parse_wcs_header(hdr))

            elif map_type == 'healpix':
                nside = hdr['NSIDE']
                nested = hdr['ORDERING'].strip().lower() in ['nest', 'nested']
                map_opts.update(nested=nested)

            # primary HDU
            if H.data is None:
                continue

            # extracting a particular HDU
            if map_type is not None and hdu is not None and hidx != hdu:
                continue

            if map_type == 'flat':
                if hdr.get('ISWEIGHT', None):
                    if 'W' not in keys:
                        continue
                else:
                    if hdr.get('POLTYPE', 'T') not in keys:
                        continue
                data = np.array(H.data, dtype=float)
                if map_opts.pop('transpose', False):
                    data = np.array(data.T)
                if uconv is not None:
                    data *= uconv

            # map type must be known if the HDU contains map data
            if map_type not in ['flat', 'healpix']:
                raise ValueError("Unknown map type in HDU {}".format(hidx))

            if map_type == 'flat' and hdr.get('ISWEIGHT', None):
                # flat map weights
                assert('T' in output)
                if 'W' not in output:
                    output['W'] = G3SkyMapWeights()
                weight_map = output['W']
                pol_type = getattr(MapPolType, hdr['WTYPE'], None)
                fm = FlatSkyMap(data, pol_type=pol_type, **map_opts)
                fm.overflow = overflow
                setattr(weight_map, hdr['WTYPE'], fm)
                del data

            elif map_type == 'flat' and not hdr.get('ISWEIGHT', None):
                # flat map data
                ptype = hdr.get('POLTYPE', 'T')
                pol_type = getattr(MapPolType, ptype, None)
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
                    'Q-POLARISATION': 'Q',
                    'U-POLARISATION': 'U',
                    'I_STOKES': 'T',
                    'Q_STOKES': 'Q',
                    'U_STOKES': 'U',
                    'I': 'T',
                    'II': 'TT',
                    'IQ': 'TQ',
                    'IU': 'TU',
                    'II_COV': 'TT',
                    'IQ_COV': 'TQ',
                    'IU_COV': 'TU',
                    'QQ_COV': 'QQ',
                    'QU_COV': 'QU',
                    'UU_COV': 'UU',
                    'I_MEAN': 'T',
                    'Q_MEAN': 'Q',
                    'U_MEAN': 'U',
                }

                pix = None

                partial = hdr.get('INDXSCHM') == 'EXPLICIT' or hdr.get('OBJECT') == 'PARTIAL'

                for cidx, hcol in enumerate(H.data.names):
                    col = col_dict.get(hcol, hcol)
                    if col in 'TQU' and col not in keys:
                        continue
                    elif col in ['TT', 'TQ', 'TU', 'QQ', 'QU', 'UU'] and 'W' not in keys:
                        continue

                    data = np.array(H.data[hcol], dtype=float).ravel()

                    if col == 'PIXEL' or (partial and cidx == 0):
                        pix = np.array(data, dtype=int).ravel()
                        del data
                        continue

                    uconv = None
                    u = 'TUNIT{:d}'.format(cidx + 1)
                    if u in hdr:
                        u = hdr[u].lower()
                        if u in unit_dict:
                            units, uconv = unit_dict[u.lower()]
                            if uconv is not None:
                                data *= uconv
                    overflow = hdr.get('TOFLW{:d}'.format(cidx + 1), overflow)

                    weighted = hdr.get(
                        'TISWGT{:d}'.format(cidx + 1),
                        col in ['TT', 'TQ', 'TU', 'QQ', 'QU', 'UU'],
                    )
                    pol_type = getattr(MapPolType, col, None)
                    mdata = (pix, data, nside) if pix is not None else data

                    if weighted:
                        assert('T' in output)
                        if 'W' not in output:
                            output['W'] = G3SkyMapWeights()
                        weight_map = output['W']
                        hm = HealpixSkyMap(mdata, pol_type=pol_type, **map_opts)
                        hm.overflow = overflow
                        setattr(weight_map, col, hm)
                    else:
                        hm = HealpixSkyMap(mdata, pol_type=pol_type, **map_opts)
                        output[col] = hm

                    del mdata, data

                del pix
                break

            del H.data

    for k, m in output.items():
        m.pol_conv = map_opts['pol_conv']
        if k == 'W':
            continue
        m.weighted = 'W' in output

    return output


def load_proj_dict(inverse=False):
    """
    Dictionary that associates projection names to their three-letter codes in WCS.
    """

    projdict = {
        'Proj0': 'SFL',
        'Proj1': 'CAR',
        'Proj2': 'SIN',
        'Proj3': 'ARC',
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
        'ProjZenithalEquidistant': 'ARC',
        'ProjARC': 'ARC',
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


def get_wcs(skymap, reset=False):
    """
    Creates WCS (world coordinate system) information from information in the
    input FlatSkyMap object.
    """

    if hasattr(skymap, '_wcs') and not reset:
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
    w.wcs.crpix = [skymap.x_center + 1, skymap.y_center + 1]
    w.wcs.crval = [
        skymap.alpha_center / core.G3Units.deg,
        skymap.delta_center / core.G3Units.deg,
    ]

    if proj_abbr in ['SFL', 'CAR', 'CEA']:
        w.wcs.crval[1] = 0.0
        if proj_abbr == 'CEA':
            v = np.sin(skymap.delta_center / core.G3Units.rad) / skymap.y_res
        else:
            v = skymap.delta_center / skymap.y_res
        w.wcs.crpix[1] -= v
        if str(skymap.proj) == 'ProjBICEP':
            w.wcs.pc[0][0] = 1. / np.cos(skymap.delta_center / core.G3Units.rad)

    skymap._wcs = w
    return w

# set property
wcs_doc = "astropy.wcs.WCS instance containing projection information"
setattr(FlatSkyMap, 'wcs', property(get_wcs, doc=wcs_doc))


def create_wcs_header(skymap):
    """
    Return a FITS header for the WCS information in the FlatSkyMap object.
    """
    header = skymap.wcs.to_header()
    header['PROJ'] = str(skymap.proj)
    header['ALPHA0'] = skymap.alpha_center / core.G3Units.deg
    header['DELTA0'] = skymap.delta_center / core.G3Units.deg
    header['X0'] = skymap.x_center
    header['Y0'] = skymap.y_center
    return header


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
    cdelt = w.wcs.get_cdelt()
    x_res = -cdelt[order[0]] * core.G3Units.deg
    res = cdelt[order[1]] * core.G3Units.deg

    # parse map center
    crval = w.wcs.crval
    alpha_center = header.get('ALPHA0', crval[order[0]]) * core.G3Units.deg
    delta_center = header.get('DELTA0', crval[order[1]]) * core.G3Units.deg

    crpix = w.wcs.crpix
    x_center = header.get('X0', crpix[order[0]] - 1)
    y_center = header.get('Y0', crpix[order[1]] - 1)

    if proj == 'ProjBICEP':
        v = 1. / np.cos(delta_center / core.G3Units.rad)
        if not np.isclose(v, w.wcs.pc[order[0]][order[0]]):
            raise ValueError("Inconsistent BICEP projection resolution")

    if wcsproj  == 'CEA':
        for i, m, v in w.wcs.get_pv():
            if i == order[1] + 1 and m == 1:
                if not np.isclose(v, 1):
                    raise ValueError(
                        "CEA projection with non-conformal scaling parameter "
                        "PV%d_1 not supported" % i
                    )

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
                     compress=False, quantize_level=16.0, hdr=None):
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
    compress : str or bool
        If defined, and if input maps are FlatSkyMap objects, store these in a
        series of compressed image HDUs, one per map.  Otherwise, store input
        maps in a series of standard ImageHDUs, which are readable with older
        FITS readers (e.g. idlastro). If defined, the compression algorithm to
        be used by the Astropy class astropy.io.fits.CompImageHDU.
        Can be: 'RICE_1', 'RICE_ONE', 'PLIO_1', 'GZIP_1', 'GZIP_2' or
        'HCOMPRESS_1'. Only GZIP_1 and GZIP_2 are lossless, although only
        for integer data.
    quantize_level : float
        Floating point quantization level for compression.  Higher values result
        in more accurate floating point representation, but worse compression
        ratio.  See the astropy FITS image documention for details:
        https://docs.astropy.org/en/stable/io/fits/api/images.html
    hdr  : dict
       If defined, extra keywords to be appened to the FITS header. The dict
       can contain entries such as ``hdr['NEWKEY'] = 'New value'`` or
       ``hdr['NEWKEY'] = ('New value', "Comment for New value")``.
    """

    import astropy.io.fits

    if isinstance(T, (FlatSkyMap, HealpixSkyMap)):
        flat = isinstance(T, FlatSkyMap)
    else:
        raise TypeError("Input map must be a FlatSkyMap or HealpixSkyMap instance")

    ctype = None
    if compress == True:
        ctype = 'GZIP_2'
    elif isinstance(compress, str):
        ctype = compress

    if Q is not None:
        assert(U is not None)
        assert(U.pol_conv == MapPolConv.IAU or U.pol_conv == MapPolConv.COSMO)
        assert(T.compatible(Q))
        assert(T.compatible(U))
        pol = True
        maps = [T, Q, U]
        if flat:
            names = ['T', 'Q', 'U']
        else:
            names = ['TEMPERATURE', 'Q_POLARISATION', 'U_POLARISATION']
    else:
        assert(U is None)
        pol = False
        maps = [T]
        names = ['T' if flat else 'TEMPERATURE']

    if flat:
        header = create_wcs_header(T)

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
    if pol:
        header['POLCCONV'] = str(U.pol_conv)
    header['POLAR'] = pol
    header['UNITS'] = str(T.units)
    header['WEIGHTED'] = T.weighted

    # Append extra metadata in hdr if defined
    if hdr:
        for keyword in hdr.keys():
            header[keyword] = hdr[keyword]

    hdulist = astropy.io.fits.HDUList()
    cols = astropy.io.fits.ColDefs([])

    try:
        for m, name in zip(maps, names):
            if flat:
                if compress:
                    hdu = astropy.io.fits.CompImageHDU(
                        np.asarray(m),
                        header=header,
                        compression_type=ctype,
                        quantize_level=quantize_level,
                    )
                else:
                    hdu = astropy.io.fits.ImageHDU(np.asarray(m), header=header)
                hdu.header['ISWEIGHT'] = False
                hdu.header['POLTYPE'] = name
                hdu.header['OVERFLOW'] = m.overflow
                hdu.header['FLATPOL'] = m.flat_pol
                hdulist.append(hdu)
                del hdu
            else:
                if not T.dense:
                    pix1, data = m.nonzero_pixels()
                    pix1 = np.asarray(pix1)
                    data = np.asarray(data)
                    data = data[np.argsort(pix1)]
                    pix1.sort()

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
                assert(W.polarized)
                wnames = ['TT', 'TQ', 'TU', 'QQ', 'QU', 'UU']
            else:
                assert(not W.polarized)
                wnames = ['TT']

            for wt in wnames:
                m = getattr(W, wt)
                assert T.compatible(m), 'Weight {} is not compatible with T'.format(wt)

                if flat:
                    if compress:
                        hdu = astropy.io.fits.CompImageHDU(
                            np.asarray(m),
                            header=header,
                            compression_type=ctype,
                            quantize_level=quantize_level,
                        )
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
                        data = data[np.argsort(pix1)]
                        pix1.sort()
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
def SaveMapFrame(
    frame,
    output_file=None,
    hdr=None,
    compress=False,
    quantize_level=16.0,
    overwrite=False,
):
    """
    Save a map frame to a FITS file.  See ``save_skymap_fits`` for details.  The
    map frame should contain T maps and (optionally) unpolarized weights, or T/Q/U
    maps and (optionally) polarized weights to store in the output file.

    Arguments
    ---------
    output_file : str or callable
        Fits filename to which the map will be written.  Maybe a callable
        function that takes a frame object as its sole input argument and
        returns a string filename.
    hdr  : dict
       If defined, extra keywords to be appened to the FITS header. The dict
       can contain entries such as ``hdr['NEWKEY'] = 'New value'`` or
       ``hdr['NEWKEY'] = ('New value', "Comment for New value")``.
    compress : str or bool
        If defined, and if input maps are FlatSkyMap objects, store these in a
        series of compressed image HDUs, one per map.  Otherwise, store input
        maps in a series of standard ImageHDUs, which are readable with older
        FITS readers (e.g. idlastro). If defined, the compression algorithm to
        be used by the Astropy class astropy.io.fits.CompImageHDU.
        Can be: 'RICE_1', 'RICE_ONE', 'PLIO_1', 'GZIP_1', 'GZIP_2' or
        'HCOMPRESS_1'. Only GZIP_1 and GZIP_2 are lossless, although only
        for integer data.
    quantize_level : float
        Floating point quantization level for compression.  Higher values result
        in more accurate floating point representation, but worse compression
        ratio.  See the astropy FITS image documention for details:
        https://docs.astropy.org/en/stable/io/fits/api/images.html
    overwrite : bool
        If True, any existing file with the same name will be ovewritten.
    """

    if frame.type != core.G3FrameType.Map:
        return

    T = frame['T']
    Q = frame.get('Q', None)
    U = frame.get('U', None)
    W = frame.get('Wpol', frame.get('Wunpol', None))

    if callable(output_file):
        output_file = output_file(frame)

    save_skymap_fits(
        output_file,
        T=T,
        Q=Q,
        U=U,
        W=W,
        overwrite=overwrite,
        compress=compress,
        quantize_level=quantize_level,
        hdr=hdr,
    )
