#!/usr/bin/env python

import numpy as np
from spt3g import core
import unittest

class TestBackAndForthConversion(unittest.TestCase):
    def test_complex_vec_conversion(self):
        c = np.arange(100000, dtype = 'complex128') + np.arange(100000, dtype = 'complex128')[::-1]*1j
        d = core.G3VectorComplexDouble(c)
        self.assertTrue(np.max( np.abs(c-d)) == 0)
        
        # Make sure item access works too; otherwise the above is essentially a test of memcpy()
        e = np.asarray(list(d))
        self.assertTrue(np.max( np.abs(c-e)) == 0)

if __name__ == '__main__':
    unittest.main()
