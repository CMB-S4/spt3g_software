#!/usr/bin/env python

from spt3g import core
import unittest

import numpy as np
import copy

SEC = core.G3Units.sec


def get_test_block(length, keys, offset=0, ordered=True):
    type_cycle = [(core.G3VectorDouble, float),
                  (core.G3VectorInt, int),
                  (core.G3VectorString, str)]
    t0 = core.G3Time('2019-01-01T12:30:00') + offset*SEC
    m = core.G3TimesampleMap()
    times = np.arange(length)
    if not ordered:
        np.random.shuffle(times)
    m.times = core.G3VectorTime([t0 + t*SEC for t in times])
    for i, k in enumerate(keys):
        y = (np.random.uniform(size=length) * 100).astype(int)
        constructor, cast_func = type_cycle[i % len(type_cycle)]
        vect = constructor(list(map(cast_func, y)))
        m[k] = vect
    return m


class TestTimesampleMap(unittest.TestCase):
    def test_00_internal_checks(self):
        # Valid block.
        m = get_test_block(100, ['x', 'y', 'z'])
        m.Check()

    def test_10_safety(self):
        m0 = get_test_block(100, ['x', 'y', 'z'])
        m1 = get_test_block(101, ['x', 'y', 'z'])
        # Try to add an incompatible element.
        with self.assertRaises(ValueError):
            m0.times = m1.times
        with self.assertRaises(ValueError):
            m0['A'] = m1['x']
        # But we should be able to change times in an empty vector.
        m0 = get_test_block(100, [])
        m1 = get_test_block(101, ['x', 'y', 'z'])
        m0.times = m1.times
        m0['x'] = m1['x']

    def test_20_concat(self):
        # Test concatenation.
        key_list = ['x', 'y', 'z']
        m0 = get_test_block(100, key_list)
        m1 = get_test_block(200, key_list, offset=100)
        m01 = m0.Concatenate(m1)
        self.assertTrue(np.all(
            np.hstack([np.array(m0.times), np.array(m1.times)]) == np.array(m01.times)))
        for k in key_list:
            self.assertTrue(np.all(
                np.hstack([np.array(m0[k]), np.array(m1[k])]) == np.array(m01[k])))
        for fail_vec in [
                get_test_block(200, key_list + ['extra'], 100),
                get_test_block(200, key_list[:-1], 100),
                ]:
            with self.assertRaises(ValueError):
                m0.Concatenate(fail_vec)

    def test_30_serialization(self):
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

    def test_40_sort(self):
        m0 = get_test_block(100, ['x', 'y', 'z'], ordered=False)
        m1 = copy.deepcopy(m0)
        m0.Sort()
        idx = np.argsort(m1.times)
        self.assertTrue((np.asarray(m1.times)[idx] == np.asarray(m0.times)).all())
        for k in ['x', 'y', 'z']:
            self.assertTrue((np.asarray(m1[k])[idx] == np.asarray(m0[k])).all())


if __name__ == '__main__':
    unittest.main()
