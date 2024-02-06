#!/usr/bin/env python

from spt3g import core

import unittest
import os

test_filename = 'inttest.g3'

# Type of int needed to store the specified value.
bit_sizes = [
    ( 8,  0),
    ( 8,  1),
    ( 8, -1),
    ( 8,  127),
    ( 8, -128),
    (16,  128),
    (16, -129),
    (16,  32767),
    (16, -32768),
    (32,  32768),
    (32, -32769),
    (32,  2147483647),
    (32, -2147483648),
    (64,  2147483648),
    (64,  2147483649),
    (64,  9223372036854775807),
    (64, -9223372036854775808),
]

class TestVectorInt(unittest.TestCase):
    def tearDown(self):
        os.remove(test_filename)
        
    def test_serialize(self):
        """Confirm full ranges can be saved and loaded."""

        w = core.G3Writer(test_filename)
        for isize, val in bit_sizes:
            f = core.G3Frame()
            f['v'] = core.G3VectorInt([val] * 10)
            w(f)
        del w

        r = core.G3Reader(test_filename)
        for isize, val in bit_sizes:
            f = r(None)[0]
            v_in = list(f['v'])
            self.assertTrue(all([_v == val for _v in v_in]),
                             "Failed to save/load value %i" % val)
        del r

    def test_serialize_map(self):
        """Confirm full ranges can be saved and loaded."""

        w = core.G3Writer(test_filename)
        for isize, val in bit_sizes:
            f = core.G3Frame()
            f['m'] = core.G3MapInt({str(i): val for i in range(10)})
            w(f)
        del w

        r = core.G3Reader(test_filename)
        for isize, val in bit_sizes:
            f = r(None)[0]
            v_in = list(f['m'].values())
            self.assertTrue(all([_v == val for _v in v_in]),
                             "Failed to save/load value %i" % val)
        del r

    def test_serialize_map_vector(self):
        """Confirm full ranges can be saved and loaded."""

        w = core.G3Writer(test_filename)
        for isize, val in bit_sizes:
            f = core.G3Frame()
            f['m'] = core.G3MapVectorInt({'a': [val] * 10})
            w(f)
        del w

        r = core.G3Reader(test_filename)
        for isize, val in bit_sizes:
            f = r(None)[0]
            v_in = list(f['m']['a'])
            self.assertTrue(all([_v == val for _v in v_in]),
                             "Failed to save/load value %i" % val)
        del r

    def test_compression(self):
        """Confirm that minimum necessary int size is used for serialization."""
        count = 10000
        overhead = 200
        for isize, val in bit_sizes:
            w = core.G3Writer(test_filename)
            f = core.G3Frame()
            f['v'] = core.G3VectorInt([val] * count)
            w(f)
            del w
            on_disk = os.path.getsize(test_filename)
            self.assertTrue(abs(on_disk - count * isize / 8) <= overhead,
                            "Storage for val %i took %.2f bytes/item, "
                            "too far from %.2f bytes/item" %
                            (val, on_disk / count, isize / 8))

    def test_compression_map(self):
        """Confirm that minimum necessary int size is used for serialization."""
        count = 10000
        overhead = 12
        for isize, val in bit_sizes:
            w = core.G3Writer(test_filename)
            f = core.G3Frame()
            f['v'] = core.G3MapInt({str(i): val for i in range(count)})
            w(f)
            del w
            on_disk = os.path.getsize(test_filename)
            self.assertTrue(abs(on_disk - count * isize / 8) <= (overhead * count),
                            "Storage for val %i took %.2f bytes/item, "
                            "too far from %.2f bytes/item" %
                            (val, on_disk / count - overhead, isize / 8))

    def test_compression_map_vector(self):
        """Confirm that minimum necessary int size is used for serialization."""
        count = 10000
        overhead = 200
        for isize, val in bit_sizes:
            w = core.G3Writer(test_filename)
            f = core.G3Frame()
            f['v'] = core.G3MapVectorInt({'a': [val] * count})
            w(f)
            del w
            on_disk = os.path.getsize(test_filename)
            self.assertTrue(abs(on_disk - count * isize / 8) <= overhead,
                            "Storage for val %i took %.2f bytes/item, "
                            "too far from %.2f bytes/item" %
                            (val, on_disk / count, isize / 8))

if __name__ == '__main__':
    unittest.main()
