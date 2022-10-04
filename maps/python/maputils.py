from spt3g import core, maps

import numpy as np

__all__ = ["flatsky_to_healpix", "healpix_to_flatsky"]


@core.usefulfunc
def healpix_to_flatsky(
    map_in, nest=False, map_stub=None, rebin=1, interp=False, **kwargs
):
    """
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
        If supplied and >1, subdivide the output pixel by n x n with each 
        sub-pixel taking on the input map values at pixel center (with interp
        or nearest neighbor). The output pixel takes on the average of the sub-pixel
        values.
        In the case that the input map has higher resolution than the output
        map (and that the input map is not low-pass filtered to remove
        information above the Nyquist freq. of the output map pixel), this
        reduces aliasing compared with direct sampling. But there would still be
        aliased power from the input map from freq above the ouput map pixel's
        Nyquist. 

    interp[false]: bool
        If True, use bilinear interpolation to extract values from the 4 closest
        pixel centers of the healpix map. Otherwise, the nearest-neighbor value 
        is used.

    All additional keyword arguments are passed to FlatSkyMap to construct the
    output map object.  Required if `map_stub` is not supplied, otherwise
    ignored.

    Maps can be rotated between Equatorial and Galactic coordinates, and/or
    change polarization convention between COSMO and IAU, by setting the
    appropriate attributes of the input and output maps.  Note that if the input
    map is a numpy array representation of a healpix map, the coordinate system
    and polarization convention are assumed to be that of the output map.
    Conversely, attributes not defined in the output map are assumed to be
    that of the input map.

    Returns:
    --------
    output_map: FlatSkyMap
        The reprojected map
    """

    # Construct the output map
    if map_stub is None:
        if len(kwargs) == 0:
            raise ValueError("Need to provide either a stub map or map parameters (projection, center, etc.)")
        if isinstance(map_in, maps.HealpixSkyMap):
            kwargs.setdefault("coord_ref", map_in.coord_ref)
            kwargs.setdefault("pol_type", map_in.pol_type)
            kwargs.setdefault("pol_conv", map_in.pol_conv)
        map_out = maps.FlatSkyMap(**kwargs)
    else:
        if not isinstance(map_stub, maps.FlatSkyMap):
            raise TypeError("Output stub must be a FlatSkyMap")
        map_out = map_stub.clone(False)

    # Populate output map pixels with interpolation and rebinning
    if not isinstance(map_in, maps.HealpixSkyMap):
        map_in = maps.HealpixSkyMap(
            map_in,
            nested=nest,
            coord_ref=map_out.coord_ref,
            weighted=map_out.weighted,
            units=map_out.units,
            pol_type=map_out.pol_type,
            pol_conv=map_out.pol_conv,
        )
    maps.reproj_map(map_in, map_out, rebin=rebin, interp=interp)

    return map_out


@core.usefulfunc
def flatsky_to_healpix(
    map_in, map_stub=None, rebin=1, interp=False, fullsky=False, **kwargs
):
    """
    Re-pixelize a map to Healpix from one of the flat projections.

    Parameters:
    -----------
    map_in: FlatSkyMap
        The input map you want to reproject

    map_stub[None]: HealpixSkyMap
        Stub output map object to be used to construct the output map.  If not
        supplied, one will be constructed using the remaining keyword arguments.

    rebin[1]: int
        If supplied and >1, subdivide the output pixel by n x n with each 
        sub-pixel taking on the input map values at pixel center (with interp
        or nearest neighbor). The output pixel takes on the average of the sub-pixel
        values.
        In the case that the input map has higher resolution than the output
        map (and that the input map is not low-pass filtered to remove
        information above the Nyquist freq. of the output map pixel), this
        reduces aliasing compared with direct sampling. But there would still be
        aliased power from the input map from freq above the ouput map pixel's
        Nyquist. 

    interp[false]: bool
        If True, use bilinear interpolation to extract values from the input
        map.  Otherwise, the nearest-neighbor value is used.

    fullsky[false]: bool
        If True a full-sky numpy array representation of the map is returned.
        Otherwise, a HealpixSkyMap instance is returned, containing only the
        pixels that overlap with the input map.

    All additional keyword arguments are passed to HealpixSkyMap to construct
    the output map object.  Required if `map_stub` is not supplied,
    otherwise ignored.

    Maps can be rotated between Equatorial and Galactic coordinates, and/or
    change polarization convention between COSMO and IAU, by setting the
    appropriate attributes of the input and output maps.  Attributes not defined
    in the output map are assumed to be that of the input map.

    Returns:
    --------
    output_map: numpy array or HealpixSkyMap
        The array containing the healpix map.  If `fullsky` is True, this is a
        numpy array, otherwise a HealpixSkyMap instance.
    """
    if not isinstance(map_in, maps.FlatSkyMap):
        raise TypeError("Input must be a FlatSkyMap")

    # Construct the output map
    if map_stub is None:
        kwargs.setdefault("coord_ref", map_in.coord_ref)
        kwargs.setdefault("pol_type", map_in.pol_type)
        kwargs.setdefault("pol_conv", map_in.pol_conv)
        map_out = maps.HealpixSkyMap(**kwargs)
    else:
        if not isinstance(map_stub, maps.HealpixSkyMap):
            raise TypeError("Output stub must be a HealpixSkyMap")
        map_out = map_stub.clone(False)

    # optimize ringsparse storage
    a0 = map_in.alpha_center
    da = map_in.x_res * map_in.shape[1]
    circ = 2 * np.pi * core.G3Units.rad
    a0 = np.mod(a0 + circ, circ)
    amin = np.mod(a0 - da + circ, circ)
    amax = np.mod(a0 + da + circ, circ)
    amin_shift = np.mod(a0 - da + circ + circ / 2, circ)
    amax_shift = np.mod(a0 + da + circ + circ / 2, circ)
    map_out.shift_ra = bool(np.abs(amax - amin) > np.abs(amax_shift - amin_shift))

    # Populate output map pixels with interpolation and rebinning
    maps.reproj_map(map_in, map_out, rebin=rebin, interp=interp)

    if fullsky:
        map_out.dense = True
        return np.asarray(map_out)
    return map_out
