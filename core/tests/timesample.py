#!/usr/bin/env python

from spt3g import core
import unittest

import numpy as np

SEC = core.G3Units.sec


def get_test_block(length, keys, offset=0):
    type_cycle = [(core.G3VectorDouble, float),
                  (core.G3VectorInt, int),
                  (core.G3VectorString, str)]
    t0 = core.G3Time('2019-01-01T12:30:00') + offset*SEC
    m = core.G3TimesampleMap()
    m.times = core.G3VectorTime([t0 + t*SEC for t in np.arange(length)])
    for i, k in enumerate(keys):
        y = (np.random.uniform(size=length) * 100).astype(int)
        constructor, cast_func = type_cycle[i % len(type_cycle)]
        vect = constructor(list(map(cast_func, y)))
        m[k] = vect
    return m


class TestIrregBlock(unittest.TestCase):
    def test_00_internal_checks(self):
        # Valid block.
        m = get_test_block(100, ['x', 'y', 'z'])
        m.Check()

        # Construct invalid blocks.
        m = core.G3TimesampleMap()
        t0 = core.G3Time('2019-01-01T12:30:00')
        m.times = core.G3VectorTime([t0, t0 + 10*SEC, t0 + 20*SEC])
        m['x'] = core.G3VectorDouble([1, 2])
        with self.assertRaises(ValueError):
            m.Check()

    def test_10_concat(self):
        # Test concatenation.
        key_list = ['x', 'y', 'z']
        m0 = get_test_block(100, key_list)
        for fail_vec in [
                get_test_block(200, key_list + ['extra'], 100),
                get_test_block(200, key_list[:-1], 100),
                ]:
            with self.assertRaises(ValueError):
                m0.Concatenate(fail_vec)

    def test_20_serialization(self):
        m0 = get_test_block(100, ['x', 'y', 'z', 'A'])
        m1 = get_test_block(200, ['x', 'y', 'z', 'A'], 100)
        m2 = m0.Concatenate(m1)
        m0.Check()
        m1.Check()
        m2.Check()
        f = core.G3Frame()
        f['irreg0'] = m0
        f['irreg1'] = m1
        core.G3Writer('test.g3').Process(f)
        f = core.G3Reader('test.g3').Process(None)[0]
        f['irreg0'].Check()
        f['irreg1'].Check()
        f['irreg0'].Concatenate(f['irreg1'])['x']


if __name__ == '__main__':
    unittest.main()
