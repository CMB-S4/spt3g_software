#!/usr/bin/env python

import unittest
import numpy as np

from spt3g.coordinateutils import (pol_angle_gal_to_eq, pol_angle_eq_to_gal,
                                   pol_qu_gal_to_eq)

from spt3g.core import G3Units

deg = G3Units.deg
rad = G3Units.rad

class TestPolarizationCoordTransform(unittest.TestCase):
    def assertAlmostEqualMod(self, val1, val2, mod=np.pi, atol=1e-7):
        self.assertTrue(
                np.isclose(val1, val2, atol=atol) or
                np.isclose(np.abs(val1-val2), mod, atol=atol/2.)
                )

    def test_galactic_to_equatorial(self):
        l_deg = 122.93191813706228
        no_rot = pol_angle_gal_to_eq(l_deg*deg/rad, 0, 0)
        self.assertAlmostEqualMod(no_rot, 0.0)
        full_rot = pol_angle_gal_to_eq(l_deg*deg/rad, np.pi/4., 0)
        self.assertAlmostEqualMod(full_rot, 0.0)
        rot_pos = pol_angle_gal_to_eq((l_deg+10.)*deg/rad, np.pi/4., 0)
        rot_neg = pol_angle_gal_to_eq((l_deg-10.)*deg/rad, np.pi/4., 0)
        self.assertAlmostEqualMod(rot_pos, -rot_neg)
        self.assertTrue(rot_neg>np.pi/2.)
        rot_opposite = pol_angle_gal_to_eq((l_deg+180.)*deg/rad, 0., 0)
        self.assertAlmostEqualMod(rot_opposite, 0)

    def test_symmetry(self):
        rot_to_eq = pol_angle_gal_to_eq(130.*deg/rad, 15.*deg/rad, 0)
        rot_to_gal = pol_angle_eq_to_gal(
                42.75267215587*deg/rad, 76.2030255724*deg/rad, rot_to_eq)
        self.assertAlmostEqualMod(rot_to_gal, 0)

        rot_to_eq = pol_angle_gal_to_eq(15.*deg/rad, 60.*deg/rad, 0)
        rot_to_gal = pol_angle_eq_to_gal(
                222.3556645733*deg/rad, 14.944650675*deg/rad, rot_to_eq)
        self.assertAlmostEqualMod(rot_to_gal, 0)

    def test_list(self):
        rot_eq_1 = pol_angle_gal_to_eq(130.*deg/rad, 15.*deg/rad, 0.)
        rot_eq_2 = pol_angle_gal_to_eq(15.*deg/rad, 60.*deg/rad, 0.)
        rot_eq_list = pol_angle_gal_to_eq([130.*deg/rad, 15.*deg/rad],
                                          [15.*deg/rad, 60.*deg/rad],
                                          [0., 0.])
        self.assertEqual(rot_eq_1, rot_eq_list[0])
        self.assertEqual(rot_eq_2, rot_eq_list[1])
        rot_gal_1 = pol_angle_eq_to_gal(130.*deg/rad, 15.*deg/rad, 0)
        rot_gal_2 = pol_angle_eq_to_gal(15.*deg/rad, 60.*deg/rad, 0)
        rot_gal_list = pol_angle_eq_to_gal([130.*deg/rad, 15.*deg/rad],
                                           [15.*deg/rad, 60.*deg/rad],
                                           [0, 0])
        self.assertEqual(rot_gal_1, rot_gal_list[0])
        self.assertEqual(rot_gal_2, rot_gal_list[1])

    def test_qu(self):
        Q_in = 1/np.sqrt(2)
        U_in = 1/np.sqrt(2)
        (Q_out, U_out) = pol_qu_gal_to_eq(130.*deg/rad, 15.*deg/rad,
                                          Q_in, U_in)
        self.assertAlmostEqual(Q_out**2+U_out**2, 1)
        self.assertNotEqual(Q_in, Q_out)

if __name__ == '__main__':
    unittest.main()
