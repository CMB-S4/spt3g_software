#!/usr/bin/env python
from spt3g import core, maps
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

ra_0, dec_0 = maps.azel.convert_azel_to_radec(az_0, el_0)
ra_1, dec_1 = maps.azel.convert_azel_to_radec(az_1, el_1)

o_az_0 = az_0 + (np.random.rand()-0.5) * 2 * core.G3Units.deg
o_el_0 = el_0 + (np.random.rand()-0.5) * 2 * core.G3Units.deg
o_ra_0, o_dec_0 = maps.azel.convert_azel_to_radec(o_az_0, o_el_0)

import astropy.coordinates
from astropy import units as u
c = astropy.coordinates.FK5(ra=np.asarray(o_ra_0)*u.rad, dec=np.asarray(o_dec_0)*u.rad)
g = c.transform_to(astropy.coordinates.Galactic())
l_test = g.l.radian
b_test = g.b.radian


def sloppy_eq(f, g, eps=1e-5):
    np.testing.assert_allclose(f, g, rtol=0, atol=eps)

def apply_trans(alpha, delta, trans):
    q_off = maps.ang_to_quat(alpha, delta)
    ra, dec = maps.quat_to_ang(trans * q_off / trans)
    if ra < 0:
        ra += 2 * np.pi
    return ra, dec

for i in range(n_samps):
    q_trans = maps.get_transform_quat(
        az_0[i], -el_0[i], ra_0[i], dec_0[i],
        az_1[i], -el_1[i], ra_1[i], dec_1[i]
    );

    t_ra_0, t_dec_0 = apply_trans(az_0[i], -el_0[i], q_trans)
    sloppy_eq([t_ra_0, t_dec_0], [ra_0[i], dec_0[i]])

    t_ra_1, t_dec_1 = apply_trans(az_1[i], -el_1[i], q_trans)
    sloppy_eq([t_ra_1, t_dec_1], [ra_1[i], dec_1[i]])

    t_ra_o, t_dec_o = apply_trans(o_az_0[i], -o_el_0[i], q_trans)
    sloppy_eq([t_ra_o, t_dec_o], [o_ra_0[i], o_dec_0[i]])

    q_trans_gal = maps.get_fk5_j2000_to_gal_quat() * q_trans
    t_l_o, t_b_o = apply_trans(o_az_0[i], -o_el_0[i], q_trans_gal)
    sloppy_eq([t_l_o, t_b_o], [l_test[i], b_test[i]])
