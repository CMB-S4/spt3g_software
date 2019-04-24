import numpy
import astropy.coordinates, astropy.units, astropy.time

from spt3g import core

spt = astropy.coordinates.EarthLocation(lat=-89.991066*astropy.units.deg,lon=-44.65*astropy.units.deg, height=2835.0*astropy.units.meter)

@core.usefulfunc
def convert_azel_to_radec(az, el, location=spt):
    '''
    When passed G3Timestreams of azimuth and elevation positions and a
    telescope location (where SPT is, by default), return an (RA, Dec)
    tuple of timestreams corresponding to the astronomical coordinates
    at which the telescope was pointing.

    Example:
    ra, dec = convert_azel_to_radec(az, el)
    '''

    assert(az.start == el.start)
    assert(az.stop == el.stop)
    assert(az.n_samples == el.n_samples)

    # record locations of bad elevation values to mark them later
    badel_inds = numpy.where((el < -90. * core.G3Units.deg) | (el > 90. * core.G3Units.deg))
    el[badel_inds] = 0. * core.G3Units.deg

    t = astropy.time.Time(numpy.asarray([i.mjd for i in az.times()]), format='mjd')
    
    k = astropy.coordinates.AltAz(az=numpy.asarray(az)/core.G3Units.deg*astropy.units.deg, alt=numpy.asarray(el)/core.G3Units.deg*astropy.units.deg, obstime=t, location=location, pressure=0)

    kt = k.transform_to(astropy.coordinates.FK5)

    ra = core.G3Timestream(numpy.asarray(kt.ra/astropy.units.deg)*core.G3Units.deg)
    dec = core.G3Timestream(numpy.asarray(kt.dec/astropy.units.deg)*core.G3Units.deg)
    dec[badel_inds] = numpy.nan

    ra.start = dec.start = az.start
    ra.stop = dec.stop = az.stop

    return (ra, dec)

@core.usefulfunc
def convert_radec_to_azel(ra, dec, location=spt):
    '''
    When passed G3Timestreams of RA and declination positions and a
    telescope location (where SPT is, by default), return an (Az, El)
    tuple of timestreams corresponding to the local coordinates
    at which the telescope was pointing.

    Example:
    az, el = convert_radec_to_azel(az, el)
    '''

    assert(ra.start == dec.start)
    assert(ra.stop == dec.stop)
    assert(ra.n_samples == dec.n_samples)

    t = astropy.time.Time(numpy.asarray([i.mjd for i in ra.times()]), format='mjd')
    
    k = astropy.coordinates.FK5(ra=numpy.asarray(ra)/core.G3Units.deg*astropy.units.deg, dec=numpy.asarray(dec)/core.G3Units.deg*astropy.units.deg)

    kt = k.transform_to(astropy.coordinates.AltAz(obstime=t, location=location, pressure=0))

    az = core.G3Timestream(numpy.asarray(kt.az/astropy.units.deg)*core.G3Units.deg)
    el = core.G3Timestream(numpy.asarray(kt.alt/astropy.units.deg)*core.G3Units.deg)

    az.start = el.start = ra.start
    az.stop = el.stop = ra.stop

    return (az, el)

@core.indexmod
def LocalToAstronomicalPointing(frame, az_timestream='BoresightAz', el_timestream='BoresightEl', ra_timestream='BoresightRa', dec_timestream='BoresightDec', Telescope=spt):
    '''
    Converts a set of timestreams in Scan frames representing Az and El pointing of the
    telescope into RA and Declination timestreams, stored in the frame under their respective
    names.
    '''

    if frame.type != core.G3FrameType.Scan:
        return

    ra, dec = convert_azel_to_radec(frame[az_timestream], frame[el_timestream], location=Telescope)

    frame[ra_timestream] = ra
    frame[dec_timestream] = dec

