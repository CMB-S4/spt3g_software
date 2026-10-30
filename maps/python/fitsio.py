from .. import core
from . import HealpixSkyMap, FlatSkyMap, G3SkyMapWeights, G3SkyMapMask
from . import MapPolType, MapPolConv, MapCoordReference, MapProjection

import numpy as np
import os
import warnings

__all__ = [
    'load_skymap_fits',
    'save_skymap_fits',
    'load_skymapmask_fits',
    'save_skymapmask_fits',
    'parse_fits_units',
    'SaveMapFrame',
]


# unit scale prefixes
unit_prefixes = {
    '': 1.0,
    'p': 1e-12,
    'n': 1e-9,
    'u': 1e-6,
    'mu': 1e-6,
    'm': 1e-3,
    'k': 1e3,
    'M': 1e6,
}

# lowercase unit suffixes, in order of matching precedence
units_types = {
    'kcmb': (core.G3TimestreamUnits.Tcmb, core.G3Units.K),
    'krj': (core.G3TimestreamUnits.Trj, core.G3Units.K),
    'k': (core.G3TimestreamUnits.Tcmb, core.G3Units.K),
    'jy/sr': (core.G3TimestreamUnits.FluxDensity, core.G3Units.Jy / core.G3Units.sr),
    'jy': (core.G3TimestreamUnits.FluxDensity, core.G3Units.Jy),
}

# approximate size in bytes of the chunks of table rows or image rows read from or
# written to disk at a time
chunk_size = 1 << 22

# map and weight component names
stokes_cols = ['T', 'Q', 'U']
weight_cols = ['TT', 'TQ', 'TU', 'QQ', 'QU', 'UU']

# healpix table column names used by other map writers
healpix_cols = {
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

@core.usefulfunc
def parse_fits_units(units):
    """
    Parse a FITS unit string into a unit type and conversion factor.

    Arguments
    ---------
    units : str
        Unit string.  Supported options are temperature units (e.g. ``'K_CMB'``,
        ``'uK_CMB'``, ``'mK_RJ'``, ``'uK'``; bare kelvin units are interpreted as CMB
        temperature units), flux density units (e.g. ``'Jy'``, ``'mJy'``,
        ``'MJy/sr'``), any of these raised to an integer power (e.g.  ``'K_CMB^2'``),
        or the name of a ``core.G3TimestreamUnits`` type (e.g.  ``'Tcmb'`` or
        ``'None'``), which denotes data already stored in native G3Units.

    Returns
    -------
    unit_type : core.G3TimestreamUnits
        Unit type
    unit_value : float
        Value of the unit in G3Units, i.e. the factor by which the input data
        should be multiplied to convert to G3Units.

    Raises
    ------
    ValueError :
        If the unit string is not recognized.
    """
    ustr = str(units).replace(' ', '').replace('_', '')

    # check if an existing timestream units type
    for uname, uval in core.G3TimestreamUnits.names.items():
        if ustr.lower() == uname.lower():
            return uval, 1.0

    # split off integer power, e.g. K_CMB^2
    base, _, power = ustr.replace('**', '^').partition('^')
    base = base.strip('()')
    try:
        power = int(power or 1)
    except ValueError:
        raise ValueError(f"Unrecognized units {units!r}") from None

    for suffix, (utype, uconv) in units_types.items():
        if base.lower().endswith(suffix):
            prefix = base[:-len(suffix)]
            if prefix in unit_prefixes:
                return utype, (unit_prefixes[prefix] * uconv) ** power

    raise ValueError(f"Unrecognized units {units!r}")


def get_map_units(hdr, ukeys, units=None, check=None):
    """
    Return the unit type and conversion factor for a map, using the first non-empty
    header entry in the list of ``ukeys``, or the parsed ``units`` tuple if the header
    units are missing or unrecognized.  If the ``check`` label is set, raise an error if
    neither is available, or if both are available and inconsistent.
    """
    hunits = next((hdr[k] for k in ukeys if str(hdr.get(k, '')).strip()), None)
    try:
        hparsed = parse_fits_units(hunits) if hunits is not None else None
    except ValueError:
        hparsed = None
        if units is not None:
            core.log_warn(
                f"Ignoring unrecognized units {hunits!r} for {check or 'map'}",
                unit='load_skymap_fits',
            )

    if check and hparsed is None and units is None:
        raise ValueError(
            f"Missing or unrecognized units {hunits!r} for {check}.  Set the units "
            "argument to define the units of the input data, e.g. units='uK_CMB'."
        )
    if check and hparsed and units:
        if (
            hparsed[0] != units[0]
            or not np.isclose(hparsed[1], units[1], rtol=1e-12, atol=0)
        ):
            raise ValueError(
                f"Units {hunits!r} for {check} are inconsistent with the units argument"
            )
    return hparsed or units


def read_table_column(filename, offset, rowdtype, nrows, name, dtype, chunk=None):
    """
    Read a column of an uncompressed FITS binary table from disk, in chunks of
    about ``chunk`` bytes (default: ``chunk_size``) of table rows, so that only
    one chunk of the table is held in memory at a time.  ``offset`` is the
    location of the table data in the file, and ``rowdtype`` is the (big-endian)
    data type of a table row.  Returns the flattened column data, converted to
    ``dtype``.
    """
    nel = int(np.prod(rowdtype[name].shape))
    out = np.empty(nrows * nel, dtype=dtype)
    step = max(1, (chunk or chunk_size) // rowdtype.itemsize)
    with open(filename, 'rb') as f:
        f.seek(offset)
        for i0 in range(0, nrows, step):
            rows = np.fromfile(f, dtype=rowdtype, count=min(step, nrows - i0))
            out[i0 * nel:(i0 + len(rows)) * nel] = rows[name].ravel()
    return out


def read_image_rows(filename, hdulist, hidx, direct=True, chunk=None):
    """
    Iterate over the image in HDU ``hidx`` of ``hdulist`` in chunks of about
    ``chunk`` bytes (default: ``chunk_size``) of rows, yielding the index of the
    first row and a float array of the rows in each chunk.  Unscaled images in
    uncompressed files (``direct``) are read from disk with plain file reads, so
    that the image is never held in memory or memory-mapped; other images are read
    through the HDU's ``section``.
    """
    import astropy.io.fits

    H = hdulist[hidx]
    nrow, ncol = H.shape
    step = max(1, (chunk or chunk_size) // (8 * ncol))
    dtypes = {8: 'u1', 16: '>i2', 32: '>i4', 64: '>i8', -32: '>f4', -64: '>f8'}
    dtype = dtypes.get(H.header['BITPIX'])
    scaled = H.header.get('BSCALE', 1) != 1 or H.header.get('BZERO', 0) != 0
    if direct and dtype and not scaled and not isinstance(H, astropy.io.fits.CompImageHDU):
        offset = hdulist.fileinfo(hidx)['datLoc']
        rowsize = ncol * np.dtype(dtype).itemsize
        with open(filename, 'rb') as f:
            for i0 in range(0, nrow, step):
                n = min(step, nrow - i0)
                f.seek(offset + i0 * rowsize)
                rows = np.fromfile(f, dtype=dtype, count=n * ncol)
                yield i0, rows.reshape(n, ncol).astype(float)
    else:
        for i0 in range(0, nrow, step):
            yield i0, np.asarray(H.section[i0:i0 + step], dtype=float)


def write_table_column(f, offset, rowdtype, name, data, ufac=1.0, chunk=None):
    """
    Write a column of a FITS binary table to the open file ``f``, in chunks of
    about ``chunk`` bytes (default: ``chunk_size``) of table rows, by reading,
    updating and rewriting each chunk of rows, so that only one chunk of the
    table is held in memory at a time.  ``offset`` is the location of the table
    data in the file, and ``rowdtype`` is the (big-endian) data type of a table
    row.  The data are divided by ``ufac`` as they are written.
    """
    step = max(1, (chunk or chunk_size) // rowdtype.itemsize)
    for i0 in range(0, len(data), step):
        n = min(step, len(data) - i0)
        f.seek(offset + i0 * rowdtype.itemsize)
        rows = np.fromfile(f, dtype=rowdtype, count=n)
        rows[name] = data[i0:i0 + n] if ufac == 1 else data[i0:i0 + n] / ufac
        f.seek(offset + i0 * rowdtype.itemsize)
        rows.tofile(f)


def is_uncompressed(filename):
    """
    Return True if the file is an uncompressed FITS file, whose data can be read
    from disk directly.
    """
    try:
        with open(filename, 'rb') as f:
            return f.read(6) == b'SIMPLE'
    except (TypeError, OSError):
        return False


def read_hdu_column(filename, hdulist, hidx, name, dtype, direct=True):
    """
    Read a column of the binary table in HDU ``hidx`` of ``hdulist``, as a
    flattened array of the given ``dtype``.  Unscaled numeric columns in
    uncompressed files (``direct``) are read from disk in chunks of rows with
    ``read_table_column``; other columns are copied out of the table by astropy.
    """
    H = hdulist[hidx]
    col = H.columns[name]
    code = str(col.format).lstrip('0123456789')[:1]
    if direct and code in 'BIJKED' and col.bscale in [None, 1] and col.bzero in [None, 0]:
        rowdtype = H.columns.dtype.newbyteorder('>')
        offset = hdulist.fileinfo(hidx)['datLoc']
        return read_table_column(
            filename, offset, rowdtype, H.header['NAXIS2'], name, dtype
        )
    return np.array(H.data[name], dtype=dtype).ravel()


def write_table(filename, cols, header, nrows, columns):
    """
    Write a FITS file containing an empty primary HDU and a binary table with the
    given column definitions, header and number of rows.  The zero-padded table
    data are written directly to disk one column at a time with
    ``write_table_column``.  ``columns`` is an iterable of (name, data, unit
    value) tuples, e.g. a generator, so that the data for each column need only
    be computed when the column is written.  A partially written file is removed
    if writing fails.
    """
    import astropy.io.fits

    hdu = astropy.io.fits.BinTableHDU.from_columns(cols, header=header, nrows=0)
    hdu.header['NAXIS2'] = nrows
    rowdtype = hdu.columns.dtype.newbyteorder('>')

    try:
        with open(filename, 'wb+') as f:
            f.write(astropy.io.fits.PrimaryHDU().header.tostring().encode('ascii'))
            f.write(hdu.header.tostring().encode('ascii'))
            offset = f.tell()
            # zero-filled table data, padded to a whole number of FITS blocks
            size = nrows * rowdtype.itemsize
            f.truncate(offset + size + (-size % 2880))
            for name, data, ufac in columns:
                write_table_column(f, offset, rowdtype, name, data, ufac)
                del data
    except BaseException:
        # don't leave a partially written file behind
        if os.path.exists(filename):
            os.remove(filename)
        raise


@core.usefulfunc
def load_skymap_fits(filename, hdu=None, keys=None, memmap=False, units=None):
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
        Argument passed to astropy.io.fits.open, which applies to maps that are
        read through astropy.  Columns of healpix tables and unscaled flat sky map
        images in uncompressed files are read directly from disk in chunks of
        rows, regardless of this setting, so that the table or image is never
        held in memory; other images are read through astropy in chunks of rows.
        Memory-mapping is not supported for scaled images.  Default: False.
    units : str, optional
        Units of the T, Q and U data stored in the file, e.g. ``'uK_CMB'`` or
        ``'K_CMB'``.  See ``parse_fits_units`` for supported options.  Required
        if the file does not contain units information for these maps (via the
        ``BUNIT`` keyword for flat sky maps, the ``TUNITn`` keywords for healpix
        maps, or the ``UNITS`` keyword written by ``save_skymap_fits``).  If the
        file contains units information that is inconsistent with this argument,
        a ValueError is raised.  Maps are converted from these units to G3Units.
        Weight maps and any other columns without units information are loaded
        unscaled.

    Returns
    -------
    a dictionary of maps keyed with e.g. 'T', 'Q', 'U' and 'W'.  The Stokes maps
    are marked as weighted according to the ``WEIGHTED`` header keyword written
    by ``save_skymap_fits``, and as unweighted if the keyword is missing.  Flat sky
    maps that were sparse when saved (``SPARSE`` header keyword) are converted
    back to sparse storage.
    """

    import astropy.io.fits

    map_type = None
    pol = None
    map_opts = {}
    output = {}

    # defaults for missing header entries
    pol_conv = None
    coord_ref = 'Equatorial'
    map_weighted = False

    if units is not None:
        units = parse_fits_units(units)

    if keys is None:
        keys = ['T', 'Q', 'U', 'W']

    # tables and images in uncompressed files can be read directly in chunks of rows
    direct = is_uncompressed(filename)

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
            if 'COORDSYS' in hdr:
                cdict = {'C': 'Equatorial', 'G': 'Galactic', 'L': 'Local'}
                coord_ref = cdict.get(hdr['COORDSYS'], coord_ref)
            else:
                coord_ref = hdr.get('COORDREF', coord_ref)
            if 'WEIGHTED' in hdr:
                map_weighted = hdr['WEIGHTED'] in [True, 'T']

            map_opts.update(
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
            if not hdr.get('NAXIS', 0):
                continue

            # extracting a particular HDU
            if map_type is not None and hdu is not None and hidx != hdu:
                continue

            # map type must be known if the HDU contains map data
            if map_type not in ['flat', 'healpix']:
                raise ValueError("Unknown map type in HDU {}".format(hidx))

            if map_type == 'flat':
                # flat map data
                ukeys = ['BUNIT', 'TUNIT', 'UNITS']
                weighted = hdr.get('ISWEIGHT', False)
                if weighted:
                    if 'W' not in keys:
                        continue
                    ptype = hdr['WTYPE']
                    wunits = core.G3TimestreamUnits.Tcmb
                    if 'T' in output:
                        wunits = output['T'].units
                    utype, uconv = get_map_units(hdr, ukeys, (wunits, 1.0))
                else:
                    ptype = hdr.get('POLTYPE', 'T')
                    if ptype not in keys:
                        continue
                    label = "{} map in HDU {} of {}".format(ptype, hidx, filename)
                    utype, uconv = get_map_units(hdr, ukeys, units, label)

                # Fill the map from chunks of image rows, transposing if necessary.
                # Maps that were sparse when saved are stored as sparse maps,
                # with only the nonzero portion of each chunk inserted into the
                # map; others are filled in place as dense maps.
                transpose = map_opts.pop('transpose', False)
                sparse = hdr.get('SPARSE', False)
                nrow, ncol = H.shape
                pol_type = getattr(MapPolType, ptype, None)
                fm = FlatSkyMap(
                    nrow if transpose else ncol,
                    ncol if transpose else nrow,
                    pol_type=pol_type,
                    units=utype,
                    **map_opts,
                )
                if not sparse:
                    fm.sparse = False
                    buf = np.asarray(fm)
                for i0, rows in read_image_rows(filename, hdulist, hidx, direct):
                    if uconv != 1:
                        rows = rows * uconv
                    # first map row and column of the chunk
                    y0, x0 = (0, i0) if transpose else (i0, 0)
                    if transpose:
                        rows = rows.T
                    if not sparse:
                        buf[y0:y0 + rows.shape[0], x0:x0 + rows.shape[1]] = rows
                        continue
                    yz = np.flatnonzero(rows.any(axis=1))
                    if not len(yz):
                        continue
                    xz = np.flatnonzero(rows.any(axis=0))
                    fm[y0 + yz[0]:y0 + yz[-1] + 1, x0 + xz[0]:x0 + xz[-1] + 1] = (
                        np.ascontiguousarray(rows[yz[0]:yz[-1] + 1, xz[0]:xz[-1] + 1])
                    )
                if sparse:
                    fm.compact()

                if weighted:
                    assert('T' in output)
                    if 'W' not in output:
                        output['W'] = G3SkyMapWeights()
                    setattr(output['W'], ptype, fm)
                else:
                    output[ptype] = fm

            elif map_type == 'healpix':
                # healpix map data
                pix = None

                partial = hdr.get('INDXSCHM') == 'EXPLICIT' or hdr.get('OBJECT') == 'PARTIAL'

                for cidx, hcol in enumerate(H.columns.names):
                    col = healpix_cols.get(hcol, hcol)
                    if col in stokes_cols and col not in keys:
                        continue
                    elif col in weight_cols and 'W' not in keys:
                        continue

                    # read one column at a time
                    if col == 'PIXEL' or (partial and cidx == 0):
                        pix = read_hdu_column(filename, hdulist, hidx, hcol, np.int64, direct)
                        continue

                    data = read_hdu_column(filename, hdulist, hidx, hcol, float, direct)

                    weighted = col in weight_cols
                    weighted = hdr.get('TISWGT{:d}'.format(cidx + 1), weighted)

                    ukeys = ['TUNIT{:d}'.format(cidx + 1), 'TUNIT', 'UNITS']
                    if col in stokes_cols:
                        label = "{} column in HDU {} of {}".format(hcol, hidx, filename)
                        utype, uconv = get_map_units(hdr, ukeys, units, label)
                    else:
                        wunits = getattr(core.G3TimestreamUnits, 'None')
                        if weighted and 'T' in output:
                            wunits = output['T'].units
                        utype, uconv = get_map_units(hdr, ukeys, (wunits, 1.0))
                    if uconv != 1:
                        data *= uconv

                    pol_type = getattr(MapPolType, col, None)
                    mdata = (pix, data, nside) if pix is not None else (data,)
                    hm = HealpixSkyMap(
                        *mdata, pol_type=pol_type, units=utype, **map_opts
                    )
                    del mdata, data

                    if weighted:
                        assert('T' in output)
                        if 'W' not in output:
                            output['W'] = G3SkyMapWeights()
                        setattr(output['W'], col, hm)
                    else:
                        output[col] = hm

                break

    for k, m in output.items():
        m.pol_conv = map_opts['pol_conv']
        if k == 'W':
            continue
        m.weighted = map_weighted

    return output


@core.usefulfunc
def load_skymapmask_fits(filename, hdu=None, memmap=False):
    """
    Load a fits file containing a binary sky map mask.  Floating pointing or integer
    data are converted to booleans.  Assumes the fits file contains a single mask,
    with flat sky data stored in a single ``ImageHDU``, and healpix masks stored in a
    single ``BinTableHDU`` in full-sky format with implicit indexing over all pixels.
    Images and tables in uncompressed files are read directly from disk in chunks of
    rows, so that only one boolean copy of the mask data is held in memory.

    Arguments
    ---------
    filename : str
        Path to fits file
    hdu : int, optional
        If supplied, the data are extract from the given HDU index.
    memmap : bool, optional
        Argument passed to astropy.io.fits.open, which applies to masks that are read
        through astropy.  Memory-mapping is not supported for scaled images.
        Default: False.

    Returns
    -------
    G3SkyMapMask :
        The binary mask.
    """
    import astropy.io.fits

    map_type = None
    coord_ref = "Equatorial"
    map_opts = {}

    # tables and images in uncompressed files can be read directly in chunks of rows
    direct = is_uncompressed(filename)

    with astropy.io.fits.open(filename, memmap=memmap) as hdulist:
        for hidx, H in enumerate(hdulist):

            hdr = H.header
            mtype = hdr.get('MAPTYPE', None)

            if mtype == 'FLAT' or 'PROJ' in hdr or 'WCSAXES' in hdr or 'CTYPE1' in hdr:
                map_type = 'flat'
            elif mtype == 'HEALPIX' or hdr.get('PIXTYPE', None) == 'HEALPIX' or 'NSIDE' in hdr:
                map_type = 'healpix'

            if 'COORDSYS' in hdr:
                cdict = {'C': 'Equatorial', 'G': 'Galactic', 'L': 'Local'}
                coord_ref = cdict.get(hdr['COORDSYS'], coord_ref)
            else:
                coord_ref = hdr.get('COORDREF', coord_ref)
            map_opts.update(coord_ref=getattr(MapCoordReference, coord_ref))

            if map_type == 'flat':
                map_opts.update(parse_wcs_header(hdr))
            elif map_type == 'healpix':
                nside = hdr['NSIDE']
                nested = hdr['ORDERING'].strip().lower() in ['nest', 'nested']
                map_opts.update(nested=nested)

            # primary HDU
            if not hdr.get('NAXIS', 0):
                continue

            # extracting a particular HDU
            if map_type is not None and hdu is not None and hidx != hdu:
                continue

            # map type must be known if the HDU contains map data
            if map_type not in ['flat', 'healpix']:
                raise ValueError("Unknown map type in HDU {}".format(hidx))

            if map_type == 'flat':
                # flat mask data, read in chunks of rows, transposing if necessary
                transpose = map_opts.pop('transpose', False)
                nrow, ncol = H.shape
                data = np.empty((ncol, nrow) if transpose else (nrow, ncol), dtype=bool)
                for i0, rows in read_image_rows(filename, hdulist, hidx, direct):
                    if transpose:
                        data[:, i0:i0 + len(rows)] = (rows != 0).T
                    else:
                        data[i0:i0 + len(rows)] = rows != 0
                parent = FlatSkyMap(data.shape[1], data.shape[0], **map_opts)

            elif map_type == 'healpix':
                # healpix mask data, stored in the first table column
                data = read_hdu_column(
                    filename, hdulist, hidx, H.columns.names[0], bool, direct
                )
                parent = HealpixSkyMap(nside, **map_opts)

            return G3SkyMapMask(parent, data.ravel())


def load_proj_dict(inverse=False):
    """
    Dictionary that associates projection names to their three-letter codes in WCS.
    """

    projdict = {
        'Proj8': '!!!',
        'ProjSFL': 'SFL',
        'ProjCAR': 'CAR',
        'ProjSIN': 'SIN',
        'ProjARC': 'ARC',
        'ProjSTG': 'STG',
        'ProjZEA': 'ZEA',
        'ProjTAN': 'TAN',
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
                     compress=False, quantize_level=16.0, hdr=None, units=None):
    """
    Save G3 map objects to a fits file.

    `FlatSkyMap` objects are stored in a series of (optionally compressed)
    `ImageHDU` entries, in which each HDU contains the projection information in
    its header in standard WCS format, along with the image data for a single
    map (one of the Stokes maps or a weight map component).  Uncompressed images
    are written to disk in chunks of rows, so that sparse maps are not converted
    to dense maps.  Compressed images are written to disk one at a time, each from
    a dense copy of the map if it is sparse, and are all stored in image extensions
    following an empty primary HDU.  Images are stored in dense format, and the
    ``SPARSE`` keyword records whether each map was sparse, so that
    ``load_skymap_fits`` can restore its storage.

    `HealpixSkyMap` objects are stored in a `BinTableHDU` extension, which
    contains the necessary header information for compatiblity with healpix map
    readers (e.g. `healpix.read_map`), and a single table with one column per
    Stokes map or weight map component.  Sparse maps are stored as cut-sky
    pixel-indexed tables, while dense maps are stored with implicit indexing
    over all pixels.  The former produces output that is equivalent to using
    `healpy.write_map` with the `partial=True` option.

    Units information is stored in the ``BUNIT`` keyword of each image HDU for
    flat sky maps, or in the ``TUNITn`` keyword of each table column for healpix
    maps.  Temperature maps are labeled ``'mK_CMB'`` or ``'mK_RJ'`` (the native
    G3Units scale) and flux density maps are converted to ``'Jy'``.  Weighted
    maps, weights and other map types are stored in native G3Units and labeled
    with the name of their ``G3TimestreamUnits`` type (e.g. ``'Tcmb'``).  Use the
    ``units`` argument to store unweighted Stokes maps in other physical units.
    The ``UNITS`` keyword stores the ``G3TimestreamUnits`` type of the T map, and
    the ``WEIGHTED`` keyword records whether the Stokes maps are weighted.  Maps
    stored without weights are treated as unweighted.

    Arguments
    ---------
    filename : str
        Path to output file.  Must not exist, unless overwrite is True.
    T, Q, U : FlatSkyMap or HealpixSkyMap
        Maps to save
    W : G3SkyMapWeights
        Weights to save with the maps.  If not supplied, the input maps are
        stored as unweighted maps.
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
    units : str, optional
        If defined, store the T, Q and U maps in these physical units (e.g.
        ``'uK_CMB'`` or ``'mJy'``) instead of the defaults described above.  See
        ``parse_fits_units`` for supported options.  The unit type must match
        that of the input maps, and the maps must be unweighted if weights are
        also supplied.
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

    # maps stored without weights are treated as unweighted
    weighted = W is not None and any(m.weighted for m in maps)
    if units is not None:
        uparsed = parse_fits_units(units)
        if weighted:
            raise ValueError(
                "Cannot store weighted maps in physical units, remove weights first"
            )

    if flat:
        header = create_wcs_header(T)

        # flat sky maps store double-precision data
        header['BITPIX'] = -64
        header['NAXIS'] = 2

        # images to write to disk: (placeholder HDU, map, unit value)
        flat_maps = []

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

        # table column definitions, and the corresponding (name, map, unit value)
        cols = []
        columns = []

        # create sorted pixel list from nonzero pixels
        pix = order = spix = None
        if T.dense:
            nrows = T.size
        else:
            pix, tdata = T.nonzero_pixels()
            pix = np.asarray(pix)
            if not np.all(pix[1:] > pix[:-1]):
                order = np.argsort(pix)
            spix = pix if order is None else pix[order]
            nrows = len(pix)

            dtype = np.min_scalar_type(-int(spix[-1] if nrows else 0) - 1)
            dtype = np.promote_types(dtype, np.int16)
            fmt = {np.dtype(np.int16): 'I', np.dtype(np.int32): 'J'}.get(dtype, 'K')
            cols.append(astropy.io.fits.Column(name='PIXEL', format=fmt))

    header['MAPTYPE'] = 'FLAT' if flat else 'HEALPIX'
    header['COORDREF'] = str(T.coord_ref)
    if pol:
        header['POLCCONV'] = str(U.pol_conv)
    header['POLAR'] = pol
    header['UNITS'] = str(T.units)
    header['WEIGHTED'] = weighted

    # Append extra metadata in hdr if defined
    if hdr:
        for keyword in hdr.keys():
            header[keyword] = hdr[keyword]

    # maps to store: (name, map, is weight)
    entries = [(name, m, False) for m, name in zip(maps, names)]
    if W is not None:
        if pol:
            assert(W.polarized)
            wnames = weight_cols
        else:
            assert(not W.polarized)
            wnames = ['TT']
        for wt in wnames:
            m = getattr(W, wt)
            assert T.compatible(m), 'Weight {} is not compatible with T'.format(wt)
            entries.append((wt, m, True))

    # units for unweighted maps
    fits_units = {
        'Tcmb': ('mK_CMB', core.G3Units.mK),
        'Trj': ('mK_RJ', core.G3Units.mK),
        'FluxDensity': ('Jy', core.G3Units.Jy),
    }

    for name, m, weight in entries:
        # weights are stored in native G3Units
        uname, ufac = str(m.units), 1.0
        if not weight:
            if units is not None:
                if uparsed[0] != m.units:
                    raise ValueError(
                        "Units {!r} incompatible with {} map units {}".format(
                            units, name, m.units
                        )
                    )
                uname, ufac = str(units).strip(), uparsed[1]
            elif uname in fits_units and not weighted:
                uname, ufac = fits_units[uname]

        if flat:
            # placeholder data, to construct the header for an image of the right
            # shape; the map data are written to disk below
            data = np.broadcast_to(np.zeros(()), m.shape)
            hdu = astropy.io.fits.ImageHDU(data, header=header)
            hdu.header['ISWEIGHT'] = weight
            hdu.header['WTYPE' if weight else 'POLTYPE'] = name
            if not weight:
                hdu.header['FLATPOL'] = m.flat_pol
            hdu.header['BUNIT'] = uname
            hdu.header['UNITS'] = str(m.units)
            hdu.header['SPARSE'] = m.sparse
            flat_maps.append((hdu, m, ufac))
            del hdu, data
        else:
            cols.append(astropy.io.fits.Column(name=name, format='D', unit=uname))
            columns.append((name, m, ufac))
            header['TISWGT{:d}'.format(len(cols))] = weight

    if overwrite:
        if os.path.exists(filename):
            os.remove(filename)
    else:
        assert(not os.path.exists(filename))

    if flat and compress:
        # Each image is compressed and written to disk one at a time, after an
        # empty primary HDU, so that at most one dense copy of a sparse map is held
        # in memory.
        try:
            astropy.io.fits.PrimaryHDU().writeto(filename)
            for hdu, m, ufac in flat_maps:
                data = np.asarray(m.clone(copy_data=True) if m.sparse else m)
                if ufac != 1:
                    data = data / ufac
                hdu = astropy.io.fits.CompImageHDU(
                    data,
                    header=hdu.header,
                    compression_type=ctype,
                    quantize_level=quantize_level,
                )
                with astropy.io.fits.open(filename, mode='append') as out:
                    out.append(hdu)
                del hdu, data
        except BaseException:
            # don't leave a partially written file behind
            if os.path.exists(filename):
                os.remove(filename)
            raise
        return

    if flat:
        # Each image is written directly to disk in chunks of rows, copied out of
        # the map, so that sparse maps are not converted to dense maps.  The first
        # image is stored in the primary HDU.
        hdulist = astropy.io.fits.HDUList()
        for hdu, _, _ in flat_maps:
            hdulist.append(hdu)
        hdulist.update_extend()
        try:
            for hdu, (_, m, ufac) in zip(hdulist, flat_maps):
                stream = astropy.io.fits.StreamingHDU(filename, hdu.header)
                step = max(1, chunk_size // (8 * m.shape[1]))
                if m.sparse:
                    # copying rows out of a sparse map has a fixed cost per copy,
                    # so copy larger blocks of rows at a time
                    step *= 4
                for i0 in range(0, m.shape[0], step):
                    if m.sparse:
                        rows = np.asarray(m[i0:i0 + step, :])
                    else:
                        rows = np.asarray(m)[i0:i0 + step]
                    stream.write(rows if ufac == 1 else rows / ufac)
                stream.close()
        except BaseException:
            # don't leave a partially written file behind
            if os.path.exists(filename):
                os.remove(filename)
            raise
        return

    # The table is written directly to disk one column at a time, in chunks of
    # rows, so that at most one temporary copy of a single map is held in memory.
    def table_columns():
        nonlocal tdata
        if pix is not None:
            yield 'PIXEL', spix, 1.0

        for idx, (name, m, ufac) in enumerate(columns):
            if pix is None:
                data = np.asarray(m)
            else:
                # make sure all of the columns have the same set of nonzero pixels
                reorder = order
                if idx == 0:
                    data, tdata = tdata, None
                else:
                    pix1, data = m.nonzero_pixels()
                    pix1 = np.asarray(pix1)
                    if not np.array_equal(pix1, pix):
                        # same pixels in a different order
                        reorder = np.argsort(pix1)
                        if not np.array_equal(pix1[reorder], spix):
                            raise ValueError(
                                "Map {} has different nonzero pixels than "
                                "the T map".format(name)
                            )
                    del pix1
                data = np.asarray(data)
                if reorder is not None:
                    data = data[reorder]
                del reorder

            yield name, data, ufac
            # free this map's data before computing the next column
            del data

    write_table(filename, cols, header, nrows, table_columns())


@core.usefulfunc
def save_skymapmask_fits(filename, mask, overwrite=False, compress=False, hdr=None):
    """
    Save G3 mask objects to a fits file.

    Masks with ``FlatSkyMap`` parents are stored in a (optionally compressed)
    ``ImageHDU`` entry, which contains the projection information in its header in
    standard WCS format, along with the image data for a single mask.  Compressed
    images follow an empty primary HDU.

    Masks with ``HealpixSkyMap`` parents are stored in a ``BinTableHDU`` extension,
    which contains the necessary header information for compatiblity with healpix map
    readers (e.g. `healpix.read_map`), and a single table with one column for the
    mask.  Masks are stored densely with implicit indexing over all pixels, and the
    table is written directly to disk in chunks of rows.

    Mask values are stored as unsigned bytes (0 or 1).

    Arguments
    ---------
    filename : str
        Path to output file.  Must not exist, unless overwrite is True.
    mask : G3SkyMapMask
        Mask to save
    overwrite : bool
        If True, any existing file with the same name will be ovewritten.
    compress : str or bool
        If defined, and if input mask has a ``FlatSkyMap`` parent, store this in a
        compressed image HDU.  Otherwise, store input mask in a standard ImageHDU,
        which is readable with older FITS readers (e.g. idlastro). If defined, the
        compression algorithm to be used by the Astropy class
        astropy.io.fits.CompImageHDU.  Can be: 'RICE_1', 'RICE_ONE', 'PLIO_1',
        'GZIP_1', 'GZIP_2' or 'HCOMPRESS_1'. Only GZIP_1 and GZIP_2 are lossless,
        although only for integer data.
    hdr  : dict
       If defined, extra keywords to be appened to the FITS header. The dict
       can contain entries such as ``hdr['NEWKEY'] = 'New value'`` or
       ``hdr['NEWKEY'] = ('New value', "Comment for New value")``.
    """
    import astropy.io.fits

    flat = isinstance(mask.parent, FlatSkyMap)

    if flat:
        header = create_wcs_header(mask.parent)
        header["BITPIX"] = 8
        header['NAXIS'] = 2

    else:
        header = astropy.io.fits.Header()

        header['PIXTYPE'] = 'HEALPIX'
        header['ORDERING'] = 'NEST' if mask.parent.nested else 'RING'
        header['NSIDE'] = mask.parent.nside
        header['INDXSCHM'] = 'IMPLICIT'
        header['OBJECT'] = 'FULLSKY'
        cdict = {'Equatorial': 'C', 'Galactic': 'G', 'Local': 'L'}
        header['COORDSYS'] = cdict[str(mask.parent.coord_ref)]
        if header['COORDSYS'] == 'C':
            header['RADESYS'] = 'FK5'
            header['EQUINOX'] = 2000.0

    header['MAPTYPE'] = 'FLAT' if flat else 'HEALPIX'
    header['COORDREF'] = str(mask.parent.coord_ref)

    # Append extra metadata in hdr if defined
    if hdr:
        for keyword in hdr.keys():
            header[keyword] = hdr[keyword]

    # masks are stored as unsigned bytes
    data = np.asarray(mask).view(np.uint8)

    if flat:
        if compress:
            ctype = None
            if compress is True:
                ctype = 'GZIP_2'
            elif isinstance(compress, str):
                ctype = compress
            # compressed images follow an empty primary HDU
            hdu = astropy.io.fits.CompImageHDU(
                data, header=header, compression_type=ctype,
            )
            hdulist = astropy.io.fits.HDUList([astropy.io.fits.PrimaryHDU(), hdu])
        else:
            # the image is stored in the primary HDU
            hdulist = astropy.io.fits.HDUList()
            hdulist.append(astropy.io.fits.ImageHDU(data, header=header))
        hdulist.writeto(filename, overwrite=overwrite)
        return

    if os.path.exists(filename):
        if not overwrite:
            raise OSError("File {!r} already exists".format(filename))
        os.remove(filename)
    col = astropy.io.fits.Column(name="MASK", format="B")
    write_table(filename, [col], header, len(data), [("MASK", data, 1.0)])


@core.indexmod
def SaveMapFrame(
    frame,
    output_file=None,
    hdr=None,
    compress=False,
    quantize_level=16.0,
    overwrite=False,
    units=None,
):
    """
    Save a map frame to a FITS file.  See ``save_skymap_fits`` for details.  The
    map frame should contain T maps and (optionally) unpolarized weights, or T/Q/U
    maps and (optionally) polarized weights to store in the output file.  Maps
    stored without weights are treated as unweighted.

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
    units : str, optional
        If defined, store the T, Q and U maps in these physical units (e.g.
        ``'uK_CMB'`` or ``'mJy'``).  The maps must be unweighted if the frame
        also contains weights.
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
        units=units,
    )
