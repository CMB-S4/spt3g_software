import numpy as np
import os
import astropy.coordinates, astropy.units, astropy.time
from astropy.utils import iers

from spt3g import core

spt = astropy.coordinates.EarthLocation(
    lat=-89.991066 * astropy.units.deg,
    lon=-44.65 * astropy.units.deg,
    height=2835.0 * astropy.units.meter,
)

try:
    if os.getenv("SPT3G_IERS_AUTO_URL"):
        iers.conf.iers_auto_url = os.getenv("SPT3G_IERS_AUTO_URL")
    if os.getenv("SPT3G_IERS_REMOTE_TIMEOUT"):
        iers.conf.remote_timeout = float(os.getenv("SPT3G_IERS_REMOTE_TIMEOUT"))
except:
    pass


@core.usefulfunc
def check_iers(g3_time):
    """
    Check whether IERS calculations will work, and load an IERS database file
    from backup if all else fails.

    Arguments
    ---------
    g3_time : G3Time instance
        Most recent time for which an IERS calculation must be computed.
    """
    t = astropy.time.Time(g3_time.mjd, format="mjd")

    # check if accessing the IERS table outright works.
    try:
        t.ut1
        return
    except:
        pass

    # if that fails, allow extrapolation
    iers.conf.auto_max_age = None
    t = astropy.time.Time(g3_time.mjd, format="mjd")
    try:
        t.ut1
        core.log_warn("IERS auto update failed, allowing extrapolation", unit="IERS")
        return
    except:
        pass

    # and if that fails, use a locally cached file that is hopefully setup correctly.
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)), "finals2000A.all")
    iers.conf.auto_download = False
    iers.IERS.iers_table = iers.IERS_A.open(fname)
    t = astropy.time.Time(g3_time.mjd, format="mjd")
    t.ut1
    core.log_warn("Using IERS table from local cache {}".format(fname), unit="IERS")


@core.usefulfunc
def convert_azel_to_radec(az, el, location=spt):
    """
    Convert timestreams of local azimuth and elevation to right ascension and
    declination.

    Arguments
    ---------
    az, el : np.ndarray or G3Timestream
        Array of local coordinates. If inputs are G3Timestream objects,
        G3Timestreams are also returned.

    Returns
    -------
    ra, dec : np.ndarray or G3Timestream
    """

    if isinstance(az, core.G3Timestream):
        assert az.start == el.start
        assert az.stop == el.stop
        assert az.n_samples == el.n_samples
    else:
        assert len(az) == len(el)

    check_iers(az.stop)

    # record locations of bad elevation values to mark them later
    badel_inds = np.where(
        (el < -90.0 * core.G3Units.deg) | (el > 90.0 * core.G3Units.deg)
    )
    el[badel_inds] = 0.0 * core.G3Units.deg

    t = astropy.time.Time(np.asarray([i.mjd for i in az.times()]), format="mjd")

    k = astropy.coordinates.AltAz(
        az=np.asarray(az) / core.G3Units.deg * astropy.units.deg,
        alt=np.asarray(el) / core.G3Units.deg * astropy.units.deg,
        obstime=t,
        location=location,
        pressure=0,
    )
    kt = k.transform_to(astropy.coordinates.FK5)

    ra = np.asarray(kt.ra / astropy.units.deg) * core.G3Units.deg
    dec = np.asarray(kt.dec / astropy.units.deg) * core.G3Units.deg
    dec[badel_inds] = np.nan

    if isinstance(az, core.G3Timestream):
        ra = core.G3Timestream(ra)
        dec = core.G3Timestream(dec)
        ra.start = dec.start = az.start
        ra.stop = dec.stop = az.stop

    return (ra, dec)


@core.usefulfunc
def convert_radec_to_azel(ra, dec, location=spt):
    """
    Convert timestreams of right ascension and declination to local
    azimuth and elevation.

    Arguments
    ---------
    ra, dec : np.ndarray or G3Timestream
        Array of Equatorial sky coordinates. If inputs are G3Timestream
        objects, G3Timestreams are also returned.

    Returns
    -------
    az, el : np.ndarray or G3Timestream
    """

    if isinstance(ra, core.G3Timestream):
        assert ra.start == dec.start
        assert ra.stop == dec.stop
        assert ra.n_samples == dec.n_samples
    else:
        assert len(ra) == len(dec)
    check_iers(ra.stop)

    t = astropy.time.Time(np.asarray([i.mjd for i in ra.times()]), format="mjd")

    k = astropy.coordinates.FK5(
        ra=np.asarray(ra) / core.G3Units.deg * astropy.units.deg,
        dec=np.asarray(dec) / core.G3Units.deg * astropy.units.deg,
    )
    kt = k.transform_to(
        astropy.coordinates.AltAz(obstime=t, location=location, pressure=0)
    )

    az = np.asarray(kt.az / astropy.units.deg) * core.G3Units.deg
    el = np.asarray(kt.alt / astropy.units.deg) * core.G3Units.deg

    if isinstance(ra, core.G3Timestream):
        az = core.G3Timestream(az)
        el = core.G3Timestream(el)
        az.start = el.start = ra.start
        az.stop = el.stop = ra.stop

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

    if isinstance(ra, core.G3Timestream):
        assert ra.start == dec.start
        assert ra.stop == dec.stop
        assert ra.n_samples == dec.n_samples
    else:
        assert len(ra) == len(dec)

    k = astropy.coordinates.FK5(
        ra=np.asarray(ra) / core.G3Units.deg * astropy.units.deg,
        dec=np.asarray(dec) / core.G3Units.deg * astropy.units.deg,
    )
    kt = k.transform_to(astropy.coordinates.Galactic)

    glon = np.asarray(kt.l / astropy.units.deg) * core.G3Units.deg
    glat = np.asarray(kt.b / astropy.units.deg) * core.G3Units.deg

    if isinstance(ra, core.G3Timestream):
        glon = core.G3Timestream(glon)
        glat = core.G3Timestream(glat)
        glon.start = glat.start = ra.start
        glon.stop = glat.stop = ra.stop

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

    if isinstance(glon, core.G3Timestream):
        assert glon.start == glat.start
        assert glon.stop == glat.stop
        assert glon.n_samples == glat.n_samples
    else:
        assert len(glon) == len(glat)

    k = astropy.coordinates.Galactic(
        l=np.asarray(glon) / core.G3Units.deg * astropy.units.deg,
        b=np.asarray(glat) / core.G3Units.deg * astropy.units.deg,
    )
    kt = k.transform_to(astropy.coordinates.FK5)

    ra = np.asarray(kt.ra / astropy.units.deg) * core.G3Units.deg
    dec = np.asarray(kt.dec / astropy.units.deg) * core.G3Units.deg

    if isinstance(ra, core.G3Timestream):
        ra = core.G3Timestream(ra)
        dec = core.G3Timestream(dec)
        ra.start = dec.start = glon.start
        ra.stop = dec.stop = glon.stop

    return (ra, dec)


@core.indexmod
def LocalToAstronomicalPointing(
    frame,
    az_timestream="BoresightAz",
    el_timestream="BoresightEl",
    ra_timestream="BoresightRa",
    dec_timestream="BoresightDec",
    Telescope=spt,
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
