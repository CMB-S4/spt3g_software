from spt3g.core import G3MapInt
import unittest

class G3MapTestSuite(unittest.TestCase):
    def setUp(self):
        self.d = {"a": 1, "b": 2, "c": 3}
        self.m = G3MapInt(self.d)

    def test_getitem(self):
        self.assertEqual(self.m["a"], 1)
        self.assertEqual(self.m["b"], 2)
        self.assertEqual(self.m["c"], 3)
        with self.assertRaises(KeyError):
            _ = self.m["d"]

    def test_setitem(self):
        self.m["a"] = 10
        self.assertEqual(self.m["a"], 10)
        self.m["d"] = 4
        self.assertEqual(self.m["d"], 4)

    def test_delitem(self):
        del self.m["b"]
        with self.assertRaises(KeyError):
            _ = self.m["b"]
        self.assertEqual(len(self.m), 2)

        del self.m["a"]
        del self.m["c"]
        self.assertTrue(not self.m)

        with self.assertRaises(KeyError):
            del self.m["a"]

    def test_contains(self):
        self.assertTrue("a" in self.m)
        self.assertTrue("b" in self.m)
        self.assertTrue("c" in self.m)
        self.assertFalse("d" in self.m)

    def test_len(self):
        self.assertEqual(len(self.m), 3)
        self.m["d"] = 4
        self.assertEqual(len(self.m), 4)
        del self.m["a"]
        self.assertEqual(len(self.m), 3)

    def test_iter(self):
        keys = sorted(list(self.m.keys()))
        expected_keys = sorted(list(self.d.keys()))
        self.assertEqual(keys, expected_keys)
        
        values = []
        for key in self.m:
            values.append(self.m[key])
        self.assertEqual(sorted(values), sorted(self.d.values()))

    def test_keys(self):
        keys_view = self.m.keys()
        self.assertEqual(sorted(list(keys_view)), sorted(list(self.d.keys())))
        self.assertEqual(len(keys_view), len(self.d))
        self.m["d"] = 4
        self.assertEqual(sorted(list(keys_view)), sorted(list(self.d.keys()) + ["d"]))

    def test_values(self):
        values_view = self.m.values()
        self.assertEqual(sorted(list(values_view)), sorted(list(self.d.values())))
        self.assertEqual(len(values_view), len(self.d))
        self.m["a"] = 10
        values = list(values_view)
        values.sort()
        expected_values = list(self.d.values())
        expected_values[0] = 10
        expected_values.sort()
        self.assertEqual(sorted(values), sorted(expected_values))

    def test_items(self):
        items_view = self.m.items()
        self.assertEqual(sorted(list(items_view)), sorted(list(self.d.items())))
        self.assertEqual(len(items_view), len(self.d))
        self.m["d"] = 4
        self.assertEqual(sorted(list(items_view)), sorted(list(self.d.items()) + [("d", 4)]))

    def test_clear(self):
        self.m.clear()
        self.assertEqual(len(self.m), 0)
        self.assertFalse(self.m)
        self.assertEqual(list(self.m.keys()), [])
        self.assertEqual(list(self.m.values()), [])
        self.assertEqual(list(self.m.items()), [])

    def test_copy(self):
        copied_dict = self.m.copy()
        self.assertEqual(len(copied_dict), len(self.m))
        self.assertEqual(sorted(list(copied_dict.items())), sorted(list(self.m.items())))
        copied_dict["a"] = 100
        self.assertNotEqual(self.m["a"], 100)
        self.assertEqual(copied_dict["a"], 100)

    def test_get(self):
        self.assertEqual(self.m.get("a"), 1)
        self.assertEqual(self.m.get("d"), None)
        self.assertEqual(self.m.get("d", 5), 5)
        self.assertEqual(self.m.get("b", 10), 2)

    def test_pop(self):
        self.assertEqual(self.m.pop("b"), 2)
        self.assertEqual(len(self.m), 2)
        with self.assertRaises(KeyError):
            _ = self.m["b"]
        self.assertEqual(self.m.pop("d", 5), 5)
        with self.assertRaises(KeyError):
            self.m.pop("d")

    def test_update(self):
        self.m.update({"b": 20, "d": 4})
        self.assertEqual(self.m["b"], 20)
        self.assertEqual(self.m["d"], 4)
        self.assertEqual(len(self.m), 4)

        self.m.update([("c", 30), ("e", 5)])
        self.assertEqual(self.m["c"], 30)
        self.assertEqual(self.m["e"], 5)
        self.assertEqual(len(self.m), 5)

        self.m.update(f=6)
        self.assertEqual(self.m["f"], 6)
        self.assertEqual(len(self.m), 6)
        
        self.m.update({})
        self.assertEqual(len(self.m), 6)

        self.m.update([])
        self.assertEqual(len(self.m), 6)
        
        self.m.update()
        self.assertEqual(len(self.m), 6)
    
    def test_bool(self):
        self.assertTrue(self.m)
        self.m.clear()
        self.assertFalse(self.m)

if __name__ == "__main__":
    unittest.main()
