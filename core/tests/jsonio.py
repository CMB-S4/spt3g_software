#!/usr/bin/env python

from spt3g import core

f = core.G3Frame(core.G3FrameType.Scan)
f["Five"] = 5
f["Fives"] = core.G3VectorInt([5, 5, 5])
f["FiveMap"] = core.G3MapInt({"five": 5})

json = f.to_json()
f2 = core.G3Frame.from_json(json)
print(json)

assert(f.type == f2.type)
for k in f:
    assert type(f[k]) == type(f2[k])
assert f["Five"] == f2["Five"]
assert (f["Fives"] == f2["Fives"]).all()
assert list(f["FiveMap"].items()) == list(f2["FiveMap"].items())
