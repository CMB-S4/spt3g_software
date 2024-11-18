import numpy as np
import os

from .. import core

__all__ = [
    "convert_azel_to_radec",
    "convert_radec_to_azel",
    "convert_radec_to_gal",
    "convert_gal_to_radec",
    "LocalToAstronomicalPointing",
    "EquatorialToGalacticPointing",
]


def get_location(location="spt"):
    """
    Return the astropy EarthLocation for the given location name.

    Arguments
    ---------
    location : str or EarthLocation instance
        If a string, must be a recognized location name.  Currently
        only "spt" is supported.
    """
    from astropy.coordinates import EarthLocation
    import astropy.units

    if isinstance(location, EarthLocation):
        return location

    location = str(location).lower()
    if location == "spt":
        return EarthLocation(
            lat=-89.991066 * astropy.units.deg,
            lon=-44.65 * astropy.units.deg,
            height=2835.0 * astropy.units.meter,
        )

    raise NotImplementedError


iers_checked = False


def check_iers(mjd):
    """
    Check whether IERS calculations will work, and load an IERS database file
    from backup if all else fails.

    Arguments
    ---------
    mjd : Array-like
        MJD timestamps for which an IERS calculation must be computed.

    Returns
    -------
    t : astropy.time.Time instance
        Time instance for the input MJD(s).
    """
    from astropy.utils import iers
    import astropy.time

    global iers_checked
    if iers_checked:
        return astropy.time.Time(mjd, format="mjd")

    try:
        if os.getenv("SPT3G_IERS_AUTO_URL"):
            iers.conf.iers_auto_url = os.getenv("SPT3G_IERS_AUTO_URL")
        if os.getenv("SPT3G_IERS_REMOTE_TIMEOUT"):
            iers.conf.remote_timeout = float(os.getenv("SPT3G_IERS_REMOTE_TIMEOUT"))
    except:
        pass

    mjd1 = np.atleast_1d(mjd)[-1]
    t1 = astropy.time.Time(mjd1, format="mjd")

    # check if accessing the IERS table outright works.
    try:
        t1.ut1
        iers_checked = True
        return astropy.time.Time(mjd, format="mjd")
    except:
        pass

    # if that fails, allow extrapolation
    iers.conf.auto_max_age = None
    t1 = astropy.time.Time(mjd1, format="mjd")
    try:
        t1.ut1
        core.log_warn("IERS auto update failed, allowing extrapolation", unit="IERS")
        iers_checked = True
        return astropy.time.Time(mjd, format="mjd")
    except:
        pass

    # and if that fails, use a locally cached file that is hopefully setup correctly.
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)), "finals2000A.all")
    iers.conf.auto_download = False
    iers.IERS.iers_table = iers.IERS_A.open(fname)
    t1 = astropy.time.Time(mjd1, format="mjd")
    t1.ut1
    core.log_warn("Using IERS table from local cache {}".format(fname), unit="IERS")
    iers_checked = True
    return astropy.time.Time(mjd, format="mjd")


def convert_deg(d, system="g3"):
    """
    Convert the input array to the appropriate system of units.

    Arguments
    ---------
    d : array_like
        Data vector
    system : str
        System to convert to, either "g3" or "astropy".  Assumes conversion from
        the other system.
    """
    import astropy.units

    system = str(system).lower()
    if system == "astropy":
        return np.asarray(d) / core.G3Units.deg * astropy.units.deg
    elif system == "g3":
        return np.asarray(d / astropy.units.deg) * core.G3Units.deg
    raise NotImplementedError


@core.usefulfunc
def convert_azel_to_radec(az, el, location="spt", mjd=None):
    """
    Convert timestreams of local azimuth and elevation to right ascension and
    declination.

    Arguments
    ---------
    az, el : np.ndarray or G3Timestream
        Array of local coordinates. If inputs are G3Timestream objects,
        G3Timestreams are also returned.
    location : str or astropy.coordinates.EarthLocation instance
        The telescope location on Earth.
    mjd : np.ndarray
        An array of times for each az/el sample.  If input az and el
        are not G3Timestreams, this argument is required.

    Returns
    -------
    ra, dec : np.ndarray or G3Timestream
    """
    import astropy.coordinates

    singleton = False
    if isinstance(az, core.G3Timestream):
        assert az.start == el.start
        assert az.stop == el.stop
        assert az.n_samples == el.n_samples
        mjd = np.asarray([i.mjd for i in az.times])
    else:
        try:
            len(az)
        except TypeError:
            singleton = True
        if singleton:
            az = np.atleast_1d(az)
            el = np.atleast_1d(el)
            mjd = np.atleast_1d(mjd)
        assert len(az) == len(el)

    t = check_iers(mjd)

    # record locations of bad elevation values to mark them later
    badel_inds = np.where(
        (el < -90.0 * core.G3Units.deg) | (el > 90.0 * core.G3Units.deg)
    )
    el[badel_inds] = 0.0 * core.G3Units.deg

    k = astropy.coordinates.AltAz(
        az=convert_deg(az, "astropy"),
        alt=convert_deg(el, "astropy"),
        obstime=t,
        location=get_location(location),
        pressure=0,
    )
    kt = k.transform_to(astropy.coordinates.FK5())

    ra = convert_deg(kt.ra, "g3")
    dec = convert_deg(kt.dec, "g3")
    dec[badel_inds] = np.nan

    if isinstance(az, core.G3Timestream):
        ra = core.G3Timestream(ra)
        dec = core.G3Timestream(dec)
        ra.start = dec.start = az.start
        ra.stop = dec.stop = az.stop
    elif singleton:
        ra = ra[0]
        dec = dec[0]

    return (ra, dec)


@core.usefulfunc
def convert_radec_to_azel(ra, dec, location="spt", mjd=None):
    """
    Convert timestreams of right ascension and declination to local
    azimuth and elevation.

    Arguments
    ---------
    ra, dec : np.ndarray or G3Timestream
        Array of Equatorial sky coordinates. If inputs are G3Timestream
        objects, G3Timestreams are also returned.
    location : str or astropy.coordinates.EarthLocation instance
        The telescope location on Earth.
    mjd : np.ndarray
        An array of times for each ra/dec sample.  If input ra and dec
        are not G3Timestreams, this argument is required.

    Returns
    -------
    az, el : np.ndarray or G3Timestream
    """
    import astropy.coordinates

    singleton = False
    if isinstance(ra, core.G3Timestream):
        assert ra.start == dec.start
        assert ra.stop == dec.stop
        assert ra.n_samples == dec.n_samples
        mjd = np.asarray([i.mjd for i in ra.times])
    else:
        try:
            len(ra)
        except TypeError:
            singleton = True
        if singleton:
            ra = np.atleast_1d(ra)
            dec = np.atleast_1d(dec)
            mjd = np.atleast_1d(mjd)
        assert len(ra) == len(dec)

    t = check_iers(mjd)

    k = astropy.coordinates.FK5(
        ra=convert_deg(ra, "astropy"), dec=convert_deg(dec, "astropy"),
    )
    kt = k.transform_to(
        astropy.coordinates.AltAz(
            obstime=t,
            location=get_location(location),
            pressure=0,
        )
    )

    az = convert_deg(kt.az, "g3")
    el = convert_deg(kt.alt, "g3")

    if isinstance(ra, core.G3Timestream):
        az = core.G3Timestream(az)
        el = core.G3Timestream(el)
        az.start = el.start = ra.start
        az.stop = el.stop = ra.stop
    elif singleton:
        az = az[0]
        el = el[0]

    return (az, el)


@core.usefulfunc
def convert_radec_to_gal(ra, dec):
    """
    Convert timestreams of right ascension and declination to Galactic
    longitude and latitude.

    Arguments
    ---------
    ra, dec : np.ndarray or G3Timestream
        Array of Equatorial sky coordinates. If inputs are G3Timestream
        objects, G3Timestreams are also returned.

    Returns
    -------
    glon, glat : np.ndarray or G3Timestream
    """
    import astropy.coordinates

    singleton = False
    if isinstance(ra, core.G3Timestream):
        assert ra.start == dec.start
        assert ra.stop == dec.stop
        assert ra.n_samples == dec.n_samples
    else:
        try:
            len(ra)
        except TypeError:
            singleton = True
        if singleton:
            ra = np.atleast_1d(ra)
            dec = np.atleast_1d(dec)
        assert len(ra) == len(dec)

    k = astropy.coordinates.FK5(
        ra=convert_deg(ra, "astropy"), dec=convert_deg(dec, "astropy"),
    )
    kt = k.transform_to(astropy.coordinates.Galactic())

    glon = convert_deg(kt.l, "g3")
    glat = convert_deg(kt.b, "g3")

    if isinstance(ra, core.G3Timestream):
        glon = core.G3Timestream(glon)
        glat = core.G3Timestream(glat)
        glon.start = glat.start = ra.start
        glon.stop = glat.stop = ra.stop
    elif singleton:
        glon = glon[0]
        glat = glat[0]

    return (glon, glat)


@core.usefulfunc
def convert_gal_to_radec(glon, glat):
    """
    Convert timestreams of Galactic longitude and latitude to right ascension
    and declination.

    Arguments
    ---------
    glon, glat : np.ndarray or G3Timestream
        Array of Galactic sky coordinates. If inputs are G3Timestream
        objects, G3Timestreams are also returned.

    Returns
    -------
    ra, dec : np.ndarray or G3Timestream
    """
    import astropy.coordinates

    singleton = False
    if isinstance(glon, core.G3Timestream):
        assert glon.start == glat.start
        assert glon.stop == glat.stop
        assert glon.n_samples == glat.n_samples
    else:
        try:
            len(glon)
        except TypeError:
            singleton = True
        if singleton:
            glon = np.atleast_1d(glon)
            glat = np.atleast_1d(glat)
        assert len(glon) == len(glat)

    k = astropy.coordinates.Galactic(
        l=convert_deg(glon, "astropy"), b=convert_deg(glat, "astropy"),
    )
    kt = k.transform_to(astropy.coordinates.FK5())

    ra = convert_deg(kt.ra, "g3")
    dec = convert_deg(kt.dec, "g3")

    if isinstance(ra, core.G3Timestream):
        ra = core.G3Timestream(ra)
        dec = core.G3Timestream(dec)
        ra.start = dec.start = glon.start
        ra.stop = dec.stop = glon.stop
    elif singleton:
        ra = ra[0]
        dec = dec[0]

    return (ra, dec)


@core.indexmod
def LocalToAstronomicalPointing(
    frame,
    az_timestream="BoresightAz",
    el_timestream="BoresightEl",
    ra_timestream="BoresightRa",
    dec_timestream="BoresightDec",
    Telescope="spt",
):
    """
    Converts a set of timestreams in Scan frames representing Az and El pointing
    of the telescope into RA and Declination timestreams, stored in the frame
    under their respective names.
    """

    if frame.type != core.G3FrameType.Scan:
        return

    ra, dec = convert_azel_to_radec(
        frame[az_timestream], frame[el_timestream], location=Telescope
    )

    frame[ra_timestream] = ra
    frame[dec_timestream] = dec


@core.indexmod
def EquatorialToGalacticPointing(
    frame,
    ra_timestream="BoresightRa",
    dec_timestream="BoresightDec",
    glon_timestream="BoresightGalLon",
    glat_timestream="BoresightGalLat",
):
    """
    Converts a set of timestreams in Scan frames representing RA and Declination
    pointing of the telescope into Galactic longitude and latitude timestreams,
    stored in the frame under their respective names.
    """

    if frame.type != core.G3FrameType.Scan:
        return

    glon, glat = convert_radec_to_gal(frame[ra_timestream], frame[dec_timestream])

    frame[glon_timestream] = glon
    frame[glat_timestream] = glat
