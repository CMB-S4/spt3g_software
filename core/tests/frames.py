#!/usr/bin/env python

from spt3g.core import G3Frame
import unittest

class TestFrames(unittest.TestCase):
    def test_put_get_delete(self):
        frame = G3Frame()
        self.assertEqual(frame.keys(), [])
        frame["key1"] = "val1"
        self.assertEqual(frame["key1"], "val1")
        del frame["key1"]
        with self.assertRaises(KeyError):
            frame["key1"]

    def test_iter(self):
        frame = G3Frame()
        self.assertEqual(frame.keys(), [])
        frame["key1"] = "val1"
        frame["key2"] = "val2"

        keys = []
        for key in frame:
            keys.append(key)
        self.assertEqual(len(keys), 2)
        self.assertTrue("key1" in keys)
        self.assertTrue("key2" in keys)

    def test_iteritems(self):
        frame = G3Frame()
        self.assertEqual(frame.keys(), [])
        frame["key1"] = "val1"
        frame["key2"] = "val2"

        keys = []
        vals = []
        for k, v in frame.iteritems():
            keys.append(k)
            vals.append(v)
        self.assertEqual(len(keys), 2)
        self.assertEqual(len(vals), 2)
        self.assertTrue("key1" in keys)
        self.assertTrue("key2" in keys)
        self.assertTrue("val1" in vals)
        self.assertTrue("val2" in vals)

    def test_values(self):
        frame = G3Frame()
        self.assertEqual(frame.keys(), [])
        frame["key1"] = "val1"
        frame["key2"] = "val2"
        keys = []
        vals = []
        self.assertTrue("val1" in frame.values())
        self.assertTrue("val2" in frame.values())
        del frame["key1"]
        del frame["key2"]
        self.assertTrue("val1" not in frame.values())

    def test_len(self):
        frame = G3Frame()
        self.assertEqual(len(frame), 0)
        frame["key1"] = "val1"
        frame["key2"] = "val2"
        self.assertEqual(len(frame), 2)
        del frame["key1"]
        del frame["key2"]
        self.assertEqual(len(frame), 0)


if __name__ == '__main__':
    unittest.main()
