#!/usr/bin/env python

import unittest
from spt3g import core
from spt3g.core import G3Units
from spt3g.coordinateutils import ChangeCoordSys
import numpy as np
import pdb

class FakeMap(core.G3Module):
    def __init__(self, test):
        super(FakeMap, self).__init__()
        self.test = test
        self.num_frames = 1
    def Process(self, frame):
        if (self.num_frames==0):
            return []
        self.num_frames -= 1
        frame = core.G3Frame(core.G3FrameType.Map)
        if self.test == 0:
            frame['ra'] = core.G3VectorDouble([0])
            frame['dec'] = core.G3VectorDouble([90*G3Units.deg])
            frame['pol'] = core.G3VectorDouble([0])
        else:
            frame['ra'] = core.G3VectorDouble([15*G3Units.deg])
            frame['dec'] = core.G3VectorDouble([45*G3Units.deg])
            frame['pol'] = core.G3VectorDouble([0])
        return frame

class TestCoordinateChange(unittest.TestCase):
    def _correct_l(self, frame):
        if frame.type == core.G3FrameType.Map:
            self.assertAlmostEqual(frame['l'][0], 192.85948098508993*G3Units.deg)

    def _correct_b(self, frame):
        if frame.type == core.G3FrameType.Map:
            self.assertAlmostEqual(frame['b'][0], 27.128251257308968*G3Units.deg)

    def test_correct_processing(self):
        pipe = core.G3Pipeline()
        pipe.Add(FakeMap, test=0)
        pipe.Add(ChangeCoordSys, coord_ref_old=core.MapCoordReference.Equatorial,
                 coord_ref_new=core.MapCoordReference.Galactic, alt_key_in="ra",
                 alt_key_out="l", az_key_in="dec", az_key_out="b", pol_key_in="pol",
                 pol_key_out="pol_gal")
        pipe.Add(self._correct_l)
        pipe.Add(self._correct_b)
        pipe.Run()

class TestCoordinatePolRotation(unittest.TestCase):
    def _correct_pol(self, frame):
        if frame.type == core.G3FrameType.Map:
            self.assertAlmostEqual(frame['pol_gal'][0], 0.03492869720186992)

    def test_correct_processing(self):
        pipe = core.G3Pipeline()
        pipe.Add(FakeMap, test=1)
        pipe.Add(ChangeCoordSys, coord_ref_old=core.MapCoordReference.Equatorial,
                 coord_ref_new=core.MapCoordReference.Galactic, alt_key_in="ra",
                 alt_key_out="l", az_key_in="dec", az_key_out="b", pol_key_in="pol",
                 pol_key_out="pol_gal")
        pipe.Add(self._correct_pol)
        pipe.Run()

if __name__ == '__main__':
    unittest.main()
