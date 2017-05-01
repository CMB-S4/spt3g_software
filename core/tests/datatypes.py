#!/usr/bin/env python

from spt3g.core import G3Bool
import unittest

class TestDatatypes(unittest.TestCase):
    def test_g3bool(self):
        t = G3Bool(True)
        f = G3Bool(False)
        self.assertTrue(t and t.value)
        self.assertFalse(f or f.value)

if __name__ == '__main__':
    unittest.main()
