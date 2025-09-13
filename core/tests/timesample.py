#!/usr/bin/env python

from spt3g import core
import unittest

import numpy as np
import copy

SEC = core.G3Units.sec



class TestTimesampleVector(unittest.TestCase):
    def test_from_list(self):
        t0 = core.G3Time('2019-01-01T12:30:00')
        vectime = core.G3VectorTime([t0, t0 + 10*SEC])
        assert(vectime[0].time == t0.time)
        assert(vectime[1].time == t0.time + 10*SEC)
    def test_from_numpy_array(self):
        t0 = core.G3Time('2019-01-01T12:30:00')
        timestamps = np.linspace(t0.time, t0.time + 1e7*SEC, 3000)
        vectime = core.G3VectorTime(timestamps)
        assert(vectime[0] == t0)
        assert(vectime[-1] == core.G3Time(timestamps[-1]))
        assert(len(vectime) == len(timestamps))
    def test_to_numpy_array(self):
        t0 = core.G3Time('2019-01-01T12:30:00')
        timestamps = np.linspace(t0.time, t0.time + 1e7*SEC, 3000, dtype='int64')
        vectime = core.G3VectorTime(timestamps)
        np.testing.assert_array_equal(vectime, timestamps)
    def test_reinit_from_numpy_array(self):
        t0 = core.G3Time('2019-01-01T12:30:00')
        timestamps = np.linspace(t0.time, t0.time + 1e7*SEC, 3000, dtype='int64')
        vectime = core.G3VectorTime(timestamps)
        revectime = core.G3VectorTime(np.asarray(vectime))
        np.testing.assert_array_equal(revectime, timestamps)
    def test_copy_constructor(self):
        t0 = core.G3VectorTime(np.array([100000000, 200000000]))
        t1 = core.G3VectorTime(t0)
        np.testing.assert_array_equal(t0, t1)


def get_test_block(length, keys=['a', 'b', 'c', 'd'],
                   offset=0, ordered=True):
    type_cycle = [(core.G3VectorDouble, float),
                  (core.G3VectorInt, int),
                  (core.G3VectorString, str),
                  (core.G3VectorBool, bool)]
    t0 = core.G3Time('2019-01-01T12:30:00') + offset*SEC
    m = core.G3TimesampleMap()
    times = np.arange(length)
    if not ordered:
        np.random.shuffle(times)
    m.times = core.G3VectorTime(t0 + times*SEC)
    for i, k in enumerate(keys):
        y = (np.random.uniform(size=length) * 100).astype(int)
        constructor, cast_func = type_cycle[i % len(type_cycle)]
        vect = constructor(list(map(cast_func, y)))
        m[k] = vect
        if not isinstance(m[k], constructor):
            raise TypeError
    return m


class TestTimesampleMap(unittest.TestCase):
    def test_00_internal_checks(self):
        # Valid block.
        m = get_test_block(100)
        m.check()

    def test_10_safety(self):
        m0 = get_test_block(100)
        m1 = get_test_block(101)
        # Try to add an incompatible element.
        with self.assertRaises(ValueError):
            m0.times = m1.times
        with self.assertRaises(ValueError):
            m0['a'] = m1['d']
        # But we should be able to change times in an empty vector.
        m0 = get_test_block(100, [])
        m1 = get_test_block(101)
        m0.times = m1.times
        m0['x'] = m1['a']

    def test_20_concat(self):
        # Test concatenation.
        key_list = ['w', 'x', 'y', 'z']
        m0 = get_test_block(100, key_list)
        m1 = get_test_block(200, key_list, offset=100)
        m01 = m0.concatenate(m1)
        self.assertTrue(np.all(
            np.hstack([np.array(m0.times), np.array(m1.times)]) == np.array(m01.times)))
        for k in key_list:
            self.assertTrue(np.all(
                np.hstack([np.array(m0[k]), np.array(m1[k])]) == np.array(m01[k])))
        for fail_vec in [
                get_test_block(200, key_list + ['extra'], 100),
                get_test_block(200, key_list[:-1], 100),
                ]:
            with self.assertRaises(RuntimeError):
                m0.concatenate(fail_vec)

    def test_30_serialization(self):
        m0 = get_test_block(100)
        m1 = get_test_block(200, offset=100)
        m2 = m0.concatenate(m1)
        m0.check()
        m1.check()
        m2.check()
        f = core.G3Frame()
        f['irreg0'] = m0
        f['irreg1'] = m1
        core.G3Writer('test.g3').Process(f)
        f = core.G3Reader('test.g3').Process(None)[0]
        f['irreg0'].check()
        f['irreg1'].check()
        f['irreg0'].concatenate(f['irreg1'])['b']

    def test_40_sort(self):
        m0 = get_test_block(100, ordered=False)
        m1 = copy.deepcopy(m0)
        m0.sort()
        idx = np.argsort(m1.times)
        self.assertTrue((np.asarray(m1.times)[idx] == np.asarray(m0.times)).all())
        for k in m0.keys():
            self.assertTrue((np.asarray(m1[k])[idx] == np.asarray(m0[k])).all())


if __name__ == '__main__':
    unittest.main()
