from spt3g import core
from spt3g.coordinateutils import HealpixSkyMap, FlatSkyMap, G3SkyMapWeights
from spt3g.coordinateutils import MapPolType, WeightType, MapCoordReference, MapProjection
from spt3g.coordinateutils import reproj_map, get_ra_dec_map_cpp

import numpy as np
import copy
import os

__all__ = [
    'get_ra_dec_map',
    'flatsky_to_healpix',
    'healpix_to_flatsky',
    'load_skymap_fits',
    'save_skymap_fits',
    'SaveMapFrame',
    'load_proj_dict',
    'create_wcs_dict',
]

def get_ra_dec_map(map_in):
    """
    Compute the position of each pixel in the map and return maps of ra and dec.

    Parameters:
    -----------
    map_in: G3SkyMap or derived object thereof
        Input map

    Returns:
    --------
    ra, dec: G3SkyMap
        Maps of the coordinates for each pixel, with the same map parameters
        as the input map.
    """

    # store in map objects
    ra = copy.copy(map_in)
    dec = copy.copy(map_in)
    get_ra_dec_map_cpp(map_in, ra, dec)

    return ra, dec


def healpix_to_flatsky(map_in, nest=False, map_stub=None, rebin=1, interp=False,
                       **kwargs):
    '''
    Re-pixelize a map from Healpix to one of the flat sky projections.

    Parameters:
    -----------
    map_in: numpy array or HealpixSkyMap
        The array containing the input healpix map to reproject.

    nest[False]: bool
        Ordering of the healpix map, if the input is a numpy array.  Ring
        ordering is assumed by default.

    map_stub[None]: FlatSkyMap
        Stub output map object to be used to construct the output map.  If not
        supplied, one will be constructed using the remaining keyword arguments.

    rebin[1]: int
        If supplied and >1, account for sub-pixel structure by integrating
        over a sub-grid on each pixel of the given dimension.  This avoids
        aliasing of power at angular scales beyond the map resolution.

    interp[false]: bool
        If True, use bilinear interpolation to extract values from the input
        map.  Otherwise, the nearest-neighbor value is used.

    **kwargs:
        All additional keyword arguments are passed to FlatSkyMap to construct
        the output map object.  Required if `map_stub` is not supplied,
        otherwise ignored.

    Returns:
    --------
    output_map: FlatSkyMap
        The reprojected map
    '''

    # Construct the output map
    if map_stub is None:
        if isinstance(map_in, HealpixSkyMap):
            kwargs.setdefault('coord_ref', map_in.coord_ref)
        map_out = FlatSkyMap(**kwargs)
    else:
        if not isinstance(map_stub, FlatSkyMap):
            raise TypeError('Output stub must be a FlatSkyMap')
        map_out = map_stub.Clone(False)

    # Populate output map pixels with interpolation and rebinning
    if not isinstance(map_in, HealpixSkyMap):
        map_in = HealpixSkyMap(map_in, nside=nside, is_nested=nest,
                               coord_ref=map_out.coord_ref)
    reproj_map(map_in, map_out, rebin=rebin, interp=interp)

    return map_out

def flatsky_to_healpix(map_in, map_stub=None, rebin=1, interp=False,
                       fullsky=False, **kwargs):
    '''
    Re-pixelize a map to Healpix from one of the flat projections.

    Parameters:
    -----------
    map_in: FlatSkyMap
        The input map you want to reproject

    map_stub[None]: HealpixSkyMap
        Stub output map object to be used to construct the output map.  If not
        supplied, one will be constructed using the remaining keyword arguments.

    rebin[1]: int
        If supplied and >1, account for sub-pixel structure by integrating
        over a sub-grid on each pixel of the given dimension.  This avoids
        aliasing of power at angular scales beyond the map resolution.

    interp[false]: bool
        If True, use bilinear interpolation to extract values from the input
        map.  Otherwise, the nearest-neighbor value is used.

    fullsky[false]: bool
        If True a full-sky numpy array representation of the map is returned.
        Otherwise, a HealpixSkyMap instance is returned, containing only the
        pixels that overlap with the input map.

    **kwargs:
        All additional keyword arguments are passed to FlatSkyMap to construct
        the output map object.  Required if `map_stub` is not supplied,
        otherwise ignored.

    Returns:
    --------
    output_map: numpy array or HealpixSkyMap
        The array containing the healpix map
        If `fullsky` is True, this is a numpy array, otherwise a
        HealpixSkyMap instance.
    '''
    if not isinstance(map_in, FlatSkyMap):
        raise TypeError('Input must be a FlatSkyMap')

    # Construct the output map
    if map_stub is None:
        kwargs.setdefault('coord_ref', map_in.coord_ref)
        map_out = HealpixSkyMap(**kwargs)
    else:
        if not isinstance(map_stub, HealpixSkyMap):
            raise TypeError('Output stub must be a HealpixSkyMap')
        map_out = map_stub.Clone(False)

    # Populate output map pixels with interpolation and rebinning
    reproj_map(map_in, map_out, rebin=rebin, interp=interp)

    if fullsky:
        map_out.dense = True
        return np.asarray(map_out)
    return map_out


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

            if mtype == 'FLAT' or 'PROJ' in hdr:
                # flat map
                if not map_type:
                    map_type = 'flat'
                elif map_type != 'flat':
                    raise ValueError(
                        "Expected a {} sky map in HDU {}, found a flat map".format(
                            map_type, hidx
                        )
                    )
            elif mtype == 'HEALPIX' or hdr.get('PIXTYPE', None) == 'HEALPIX':
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
                proj = hdr.get('PROJ', proj)
                alpha_center = hdr.get('ALPHA0', alpha_center)
                delta_center = hdr.get('DELTA0', delta_center)
                res = hdr.get('RES', res)
                x_res = hdr.get('XRES', hdr.get('RES', xres))

                map_opts.update(
                    proj=getattr(MapProjection, proj),
                    alpha_center=alpha_center,
                    delta_center=delta_center,
                    res=res,
                    x_res=x_res,
                )

            elif map_type == 'healpix':
                nside = hdr['NSIDE']
                map_opts.update(is_nested=(hdr['ORDERING'] == 'nest'))

            # primary HDU
            if H.data is None:
                continue

            # extracting a particular HDU
            if map_type is not None and hdu is not None and hidx != hdu:
                continue

            data = np.array(H.data)

            # map type must be known if the HDU contains map data
            if map_type not in ['flat', 'healpix']:
                raise ValueError("Unknown map type in HDU {}".format(hidx))

            if map_type == 'flat' and hdr['ISWEIGHT']:
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

            elif map_type == 'flat' and not hdr['ISWEIGHT']:
                # flat map data
                pol_type = getattr(MapPolType, hdr['POLTYPE'])
                if pol_type == MapPolType.U and polcconv == 'COSMO':
                    data *= -1
                fm = FlatSkyMap(data, pol_type=pol_type, **map_opts)
                fm.overflow = overflow
                output[hdr['POLTYPE']] = fm

            elif map_type == 'healpix':
                # healpix map data

                col_dict = {
                    'I_STOKES': 'T',
                    'Q_STOKES': 'Q',
                    'U_STOKES': 'U',
                    'I': 'T',
                    'II': 'TT',
                    'IQ': 'TQ',
                    'IU': 'TU',
                }

                pix = None

                partial = hdr.get('INDXSCHM') == 'EXPLICIT' or hdr.get('OBJECT') == 'PARTIAL'

                for cidx, hcol in enumerate(H.data.names):
                    col = col_dict.get(hcol, hcol)
                    data = np.array(H.data[hcol], dtype=float)

                    if col == 'PIXEL' or (partial and cidx == 0):
                        pix = np.array(data, dtype=int)
                        continue

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

                        if pix is not None:
                            data = (pix, data, nside)
                        hm = HealpixSkyMap(data, **map_opts)
                        hm.overflow = overflow

                        setattr(weight_map, col, hm)

                    else:
                        pol_type = getattr(MapPolType, col)
                        if pol_type == MapPolType.U and polcconv == 'COSMO':
                            data *= -1
                        if pix is not None:
                            data = (pix, data, nside)

                        hm = HealpixSkyMap(data, pol_type=pol_type, **map_opts)
                        output[col] = hm

                break

    for k, m in output.items():
        if k == 'W':
            continue
        m.is_weighted = 'W' in output

    return output


@core.usefulfunc
def load_proj_dict():
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
        'ProjHealpix': '!!!',
    }

    return projdict


@core.usefulfunc
def create_wcs_dict(skymap):
    """
    Creates WCS (world coordinate system) information from information in the
    input G3SkyMap object.
    """

    projdict = load_proj_dict()
    proj_abbr = projdict[str(skymap.proj)]

    wcsdict = {
        'CRVAL1': skymap.alpha_center / core.G3Units.deg,
        'CRVAL2': skymap.delta_center / core.G3Units.deg,
        'CRPIX1': skymap.shape[1] / 2.0 + 0.5,
        'CRPIX2': skymap.shape[0] / 2.0 + 0.5,
        'CD1_1': -skymap.x_res / core.G3Units.deg,
        'CD2_2': skymap.res / core.G3Units.deg,
        'CD2_1': 0.,
        'CD1_2': 0.,
        'CTYPE1': 'RA---' + proj_abbr,
        'CTYPE2': 'DEC--' + proj_abbr,
        'CUNIT1': 'DEG',
        'CUNIT2': 'DEG',
    }

    # special voodoo for proj0 and proj1
    if str(skymap.proj) in ['Proj0', 'Proj1', 'ProjSansonFlamsteed']:
        wcsdict['CRVAL2'] = 0.
        wcsdict['CRPIX2'] -= skymap.delta_center / skymap.res

    return wcsdict


@core.usefulfunc
def save_skymap_fits(filename, T, Q=None, U=None, W=None, overwrite=False,
                     primary_header=True, compress=True):
    """
    Save G3 map objects to a fits file.

    Arguments
    ---------
    filename : str
        Path to output file.  Must not exist, unless overwrite is True.
    T[, Q, U] : G3SkyMap
        Maps to save
    W : G3SkyMapWeights
        Weights to save with the maps
    overwrite : bool
        If True, any existing file with the same name will be ovewritten.
    primary_header : bool
        If True, use the first (primary) HDU as a header, and store
        map data in extensions.  Otherwise, the T map is written directly to the
        primary HDU.  Use False to store a T-only map in an IDL-compatible file format.
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

    header = astropy.io.fits.Header()
    header['MAPTYPE'] = 'FLAT' if flat else 'HEALPIX'

    header['COORDREF'] = str(T.coord_ref)
    header['RADESYS'] = 'FK5'
    header['EQUINOX'] = '2000'
    header['POLCCONV'] = 'IAU'
    header['UNITS'] = str(T.units)
    header['WEIGHTED'] = T.is_weighted
    if W is not None:
        header['WGTTYPE'] = str(W.weight_type)

    if flat:
        header['XDIM'] = T.shape[1]
        header['YDIM'] = T.shape[0]
        header['PROJ'] = str(T.proj)
        header['ALPHA0'] = T.alpha_center
        header['DELTA0'] = T.delta_center
        header['RES'] = T.res
        header['XRES'] = T.x_res

        wcs = create_wcs_dict(T)
        for k, v in wcs.items():
            header[k] = v

        bitpix = {
            np.dtype(np.int16): 16,
            np.dtype(np.int32): 32,
            np.dtype(np.float32): -32,
            np.dtype(np.float64): -64,
        }
        header['BITPIX'] = bitpix[np.asarray(T).dtype]
        header['NAXIS'] = 2

    else:
        header['PIXTYPE'] = 'HEALPIX'
        header['ORDERING'] = 'NEST' if T.nested else 'RING'
        header['NSIDE'] = T.nside
        if T.dense:
            header['INDXSCHM'] = 'IMPLICIT'
            header['OBJECT'] = 'FULLSKY'
        else:
            header['INDXSCHM'] = 'EXPLICIT'
            header['OBJECT'] = 'PARTIAL'

        conv = {
            np.dtype(np.int16): 'I',
            np.dtype(np.int32): 'J',
            np.dtype(np.int64): 'K',
            np.dtype(np.float32): 'E',
            np.dtype(np.float64): 'D',
        }
        cols = []
        pix = None

    hdulist = astropy.io.fits.HDUList()

    if primary_header:
        hdu = astropy.io.fits.PrimaryHDU(header=header)
        hdulist.append(hdu)

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
                    cols.append(col)
                else:
                    # XXX different sparse maps can have different nonzero pixels
                    assert(len(pix) == len(pix1) and (pix == pix1).all())

            else:
                data = np.asarray(m)

            fmt = conv.get(data.dtype, 'D')
            col = astropy.io.fits.Column(
                name=name, format=fmt, array=data, unit=str(m.units)
            )
            cols.append(col)

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
            else:
                if not T.dense:
                    pix1, data = m.nonzero_pixels()
                    data = np.asarray(data)
                    assert(len(pix) == len(pix1) and (pix == pix1).all())
                else:
                    data = np.asarray(m)

                fmt = conv.get(data.dtype, 'D')
                col = astropy.io.fits.Column(
                    name=wt, format=fmt, array=data, unit=str(m.units)
                )
                cols.append(col)

                idx = len(cols)
                header['TISWGT{:d}'.format(idx)] = False
                header['TOFLW{:d}'.format(idx)] = m.overflow

    if not flat:
        hdu = astropy.io.fits.BinTableHDU.from_columns(cols, header=header)
        hdulist.append(hdu)

    if overwrite:
        if os.path.exists(filename):
            os.remove(filename)
    else:
        assert(not os.path.exists(filename))

    hdulist.writeto(filename)


@core.indexmod
def SaveMapFrame(frame, map_id, output_file, overwrite=False, primary_header=True):
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
        primary_header=primary_header,
    )
