from spt3g import core

from spt3g.coordinateutils import CutSkyHealpixMap, HealpixHitPix, FlatSkyMap
from spt3g.coordinateutils import reproj_map, reproj_fullsky_healpix_map, get_ra_dec_map_cpp

import numpy as np
import copy

__all__ = ['get_ra_dec_map', 'flatsky_to_healpix', 'healpix_to_flatsky']

def get_ra_dec_map(map_in):
    """
    Compute the position of each pixel in the map and return maps of ra and dec.

    Parameters:
    -----------
    map_in: spt3g.core.G3SkyMap or derived object thereof
        Input map

    Returns:
    --------
    ra, dec: spt3g.core.G3SkyMap
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
    map_in: numpy array or coordinateutils.CutSkyHealpixMap
        The array containing the input healpix map to reproject.

    nest[False]: bool
        Ordering of the healpix map, if the input is a numpy array.  Ring
        ordering is assumed by default.

    map_stub[None]: coordinateutils.FlatSkyMap
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
        All additional keyword arguments are passed to
        coordinateutils.FlatSkyMap to construct the output map object.  Required
        if `map_stub` is not supplied, otherwise ignored.

    Returns:
    --------
    output_map: coordinateutils.FlatSkyMap
        The reprojected map
    '''

    # Construct the output map
    if map_stub is None:
        if isinstance(map_in, CutSkyHealpixMap):
            kwargs.setdefault('coord_ref', map_in.coord_ref)
        map_out = FlatSkyMap(**kwargs)
    else:
        if not isinstance(map_stub, FlatSkyMap):
            raise TypeError('Output stub must be a FlatSkyMap')
        map_out = copy.copy(map_stub)

    # Populate output map pixels with interpolation and rebinning
    if isinstance(map_in, CutSkyHealpixMap):
        reproj_map(map_in, map_out, rebin=rebin, interp=interp)
    else:
        reproj_fullsky_healpix_map(
            map_in, map_out, nest=nest, rebin=rebin, interp=interp
        )

    return map_out

def flatsky_to_healpix(map_in, nside, nest=False, rebin=1, interp=False,
                       fullsky=False):
    '''
    Re-pixelize a map to Healpix from one of the flat projections.

    Parameters:
    -----------
    map_in: coordinateutils.FlatSkyMap
        The input map you want to reproject

    nside: int
        The nside of the output healpix map.
        nsides larger than 2048 not recommended.

    nest[False]: bool
        Ordering of the healpix map. 'Ring' by default

    rebin[1]: int
        If supplied and >1, account for sub-pixel structure by integrating
        over a sub-grid on each pixel of the given dimension.  This avoids
        aliasing of power at angular scales beyond the map resolution.

    interp[false]: bool
        If True, use bilinear interpolation to extract values from the input
        map.  Otherwise, the nearest-neighbor value is used.

    fullsky[false]: bool
        If True a full-sky numpy array representation of the map is returned.
        Otherwise, a CutSkyHealpixMap instance is returned, containing only the
        pixels that overlap with the input map.

    Returns:
    --------
    output_map: numpy array or CutSkyHealpixMap
        The array containing the healpix map
        If `fullsky` is True, this is a numpy array, otherwise a
        CutSkyHealpixMap instance.
    '''
    if not isinstance(map_in, FlatSkyMap):
        raise TypeError('Input must be a FlatSkyMap')

    # Construct the output map
    hitpix = HealpixHitPix(map_in, nside=nside, is_nested=nest)
    map_out = CutSkyHealpixMap(hitpix)

    # Populate output map pixels with interpolation
    reproj_map(map_in, map_out, rebin=rebin, interp=interp)

    if fullsky:
        return np.asarray(map_out.get_fullsky_map())
    return map_out
