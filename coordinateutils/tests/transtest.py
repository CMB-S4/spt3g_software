#!/usr/bin/env python
from spt3g import core, coordinateutils
import numpy as np

np.random.seed(42)

n_samps = int(1e3)

az_0 = core.G3Timestream(np.random.rand(n_samps)*2.0*np.pi - np.pi)
pole_avoidance = 0.7
el_0 = core.G3Timestream(np.random.rand(n_samps)*np.pi*pole_avoidance-np.pi/2.0*pole_avoidance)

az_0.start = core.G3Time('20170329_000001')
el_0.start = core.G3Time('20170329_000001')

az_0.stop = core.G3Time('20170329_100001')
el_0.stop = core.G3Time('20170329_100001')

az_1 = az_0 + 10 * core.G3Units.arcmin
el_1 = el_0 + 10 * core.G3Units.arcmin

ra_0, dec_0 = coordinateutils.azel.convert_azel_to_radec(az_0, el_0)
ra_1, dec_1 = coordinateutils.azel.convert_azel_to_radec(az_1, el_1)

o_az_0 = az_0 + (np.random.rand()-0.5) * 2 * core.G3Units.deg
o_el_0 = el_0 + (np.random.rand()-0.5) * 2 * core.G3Units.deg
o_ra_0, o_dec_0 = coordinateutils.azel.convert_azel_to_radec(o_az_0, o_el_0)

import astropy.coordinates
from astropy import units as u
c = astropy.coordinates.FK5(ra=np.asarray(o_ra_0)*u.rad, dec=np.asarray(o_dec_0)*u.rad)
g = c.transform_to(astropy.coordinates.Galactic)
l_test = g.l.radian
b_test = g.b.radian


def wrap_ring(a):
    np.asarray(a)[np.where(a > np.pi)] -= 2*np.pi
def sloppy_eq(f,g, eps = 1e-5):
    return np.abs(f-g) < eps


for i in range(n_samps):
    t_ra_0, t_dec_0 = coordinateutils.test_trans( az_0[i], el_0[i], ra_0[i], dec_0[i], 
                                            az_1[i], el_1[i], ra_1[i], dec_1[i],
                                            az_0[i], el_0[i])
    if t_ra_0 <0:
        t_ra_0 += 2*np.pi
    if not sloppy_eq(t_ra_0, ra_0[i]):
        print(t_ra_0, ra_0[i])
        print(az_0[i], el_0[i], ra_0[i], dec_0[i],
              az_1[i], el_1[i], ra_1[i], dec_1[i],
              az_0[i], el_0[i])
        assert(0)
    if not sloppy_eq(t_dec_0, dec_0[i]):
        print(t_dec_0, dec_0[i])
        print(az_0[i], el_0[i], ra_0[i], dec_0[i],
              az_1[i], el_1[i], ra_1[i], dec_1[i],
              az_0[i], el_0[i])
        assert(0)        

    t_ra_1, t_dec_1 = coordinateutils.test_trans( az_0[i], el_0[i], ra_0[i], dec_0[i], 
                                                   az_1[i], el_1[i], ra_1[i], dec_1[i],
                                                   az_1[i], el_1[i])
    if t_ra_1 <0:
        t_ra_1 += 2*np.pi

    if not sloppy_eq(t_ra_1, ra_1[i]):
        print(t_ra_1, ra_1[i])
        print(az_0[i], el_0[i], ra_0[i], dec_0[i],
              az_1[i], el_1[i], ra_1[i], dec_1[i])
        assert(0)
    if not sloppy_eq(t_dec_1, dec_1[i]):
        print(t_dec_1, dec_1[i])
        print(az_0[i], el_0[i], ra_0[i], dec_0[i],
              az_1[i], el_1[i], ra_1[i], dec_1[i])
        assert(0)        

    t_ra_o, t_dec_o = coordinateutils.test_trans( az_0[i], el_0[i], ra_0[i], dec_0[i], 
                                           az_1[i], el_1[i], ra_1[i], dec_1[i],
                                           o_az_0[i], o_el_0[i])
    if t_ra_o <0:
        t_ra_o += 2*np.pi
    if not sloppy_eq(t_ra_o, o_ra_0[i]):
        print(t_ra_o, o_ra_0[i])
        assert(0)
    if not sloppy_eq(t_dec_o, o_dec_0[i]):
        print(t_dec_o, o_dec_0[i])
        assert(0)        


        
    t_l_o, t_b_o = coordinateutils.test_gal_trans( az_0[i], el_0[i], ra_0[i], dec_0[i], 
                                            az_1[i], el_1[i], ra_1[i], dec_1[i],
                                            o_az_0[i], o_el_0[i])

    if t_l_o < 0:
        t_l_o += 2*np.pi
    if not sloppy_eq(t_l_o, l_test[i]):
        assert(0)
    if not sloppy_eq(t_b_o, b_test[i]):
        assert(0)

