#!/usr/bin/env python

from spt3g.core import G3FrameType
import unittest

class TestFramesTypes(unittest.TestCase):
    def test_key(self):
        self.assertEqual(G3FrameType.Map.key, 'M')

if __name__ == '__main__':
    unittest.main()
