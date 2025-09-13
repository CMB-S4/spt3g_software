from spt3g import core
import numpy as np

vi = core.G3VectorInt([1, 2, 3])
mi = core.G3MapInt({"a": 1, "b": 2, "c": 3})

# test default constructor
m = core.G3MapFrameObject()

# test modifiers
m["v"] = vi
assert isinstance(m["v"], core.G3VectorInt)
np.testing.assert_array_equal(m["v"], vi)
m.update(m=mi)
assert isinstance(m["m"], core.G3MapInt)
np.testing.assert_array_equal(list(m["m"].keys()), list(mi.keys()))
np.testing.assert_array_equal(list(m["m"].values()), list(mi.values()))

# test constructor from iterable
m2 = core.G3MapFrameObject({"v": vi, "m": mi})
assert isinstance(m2["v"], core.G3VectorInt)
assert isinstance(m2["m"], core.G3MapInt)
np.testing.assert_array_equal(list(m.keys()), list(m2.keys()))
np.testing.assert_array_equal(list(m.values()), list(m2.values()))

# test copy constructor
m3 = core.G3MapFrameObject(m2)
assert isinstance(m3["v"], core.G3VectorInt)
assert isinstance(m3["m"], core.G3MapInt)
np.testing.assert_array_equal(list(m3.keys()), list(m2.keys()))
np.testing.assert_array_equal(list(m3.values()), list(m2.values()))

# test default constructor
v = core.G3VectorFrameObject()

# test modifiers
v.append(vi)
assert isinstance(v[0], core.G3VectorInt)
np.testing.assert_array_equal(v[0], vi)
v.insert(0, mi)
assert isinstance(v[0], core.G3MapInt)
np.testing.assert_array_equal(list(v[0].keys()), list(mi.keys()))
np.testing.assert_array_equal(list(v[0].values()), list(mi.values()))
v[0] = vi
assert isinstance(v[0], core.G3VectorInt)
v[1] = mi
assert isinstance(v[1], core.G3MapInt)

# test constructor from iterable
v2 = core.G3VectorFrameObject([vi, mi])
assert isinstance(v2[0], core.G3VectorInt)
assert isinstance(v2[1], core.G3MapInt)
np.testing.assert_array_equal(v2, v)

# test copy constructor
v3 = core.G3VectorFrameObject(v2)
assert isinstance(v3[0], core.G3VectorInt)
assert isinstance(v3[1], core.G3MapInt)
np.testing.assert_array_equal(v3, v)
