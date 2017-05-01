#!/usr/bin/env python

from spt3g.core import G3FrameType
import unittest

class TestFramesTypes(unittest.TestCase):
    def test_from_str(self):
        types = G3FrameType.from_string("HO")
        self.assertTrue(G3FrameType.Housekeeping in types)
        self.assertTrue(G3FrameType.Observation in types)
        self.assertTrue(G3FrameType.Map not in types)

    def test_key(self):
        self.assertEqual(G3FrameType.Map.key, 'M')

if __name__ == '__main__':
    unittest.main()
