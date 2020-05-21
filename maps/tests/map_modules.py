from spt3g import core, maps
import numpy as np

pipe = core.G3Pipeline()

m = maps.FlatSkyMap(300, 300, core.G3Units.arcmin, proj=maps.MapProjection.ProjZEA)

pipe.Add(core.G3InfiniteSource, type=core.G3FrameType.Observation, n=1)
pipe.Add(maps.InjectMapStub, map_id="test_map", map_stub=m, polarized=False, weighted=True)

def RandomMap(frame):
    if frame.type != core.G3FrameType.Map:
        return

    tmap = frame.pop("T")
    tmap[:] = np.random.randn(*m.shape)
    frame["T"] = tmap

    wmap = frame.pop("Wunpol")
    wmap.TT[:] = 10 * np.ones(m.shape)
    frame["Wunpol"] = wmap


pipe.Add(RandomMap)

mex0 = maps.ExtractMaps(copy=True)
pipe.Add(mex0)

pipe.Add(maps.MakeMapsPolarized)
pipe.Add(maps.RemoveWeights)
pipe.Add(maps.ApplyWeights)

mex1 = maps.ExtractMaps()
pipe.Add(mex1)

pipe.Add(
    maps.ReplicateMaps,
    input_map_id="test_map",
    output_map_ids=["test_map1", "test_map2"],
    copy_weights=True,
)

mex2 = maps.ExtractMaps()
pipe.Add(mex2)

pipe.Add(maps.MakeMapsUnpolarized)

mex3 = maps.ExtractMaps()
pipe.Add(mex3)

tmap = m.clone(False)
tmap[:] = np.random.randn(*m.shape)
tmap.pol_type = maps.MapPolType.T
tmap.weighted = False
pipe.Add(maps.InjectMaps, map_id="test_map", maps_in=[tmap])

pipe.Run()

# check that mex0 and mex1 are nearly (but not exactly) identical
assert np.max(mex0.maps["test_map"]["T"] - mex1.maps["test_map"]["T"]) != 0
assert np.allclose(mex0.maps["test_map"]["T"], mex1.maps["test_map"]["T"])

# check that replication dropped the test map
assert "test_map" not in mex2.maps
assert "test_map" not in mex3.maps

# check that replicated maps are properly populated
for map_id in ["test_map1", "test_map2"]:
    assert map_id in mex2.maps
    assert map_id in mex3.maps

    mdict = mex2.maps[map_id]
    for k in ["T", "Q", "U", "Wpol"]:
        assert k in mdict
        assert mdict[k].compatible(m)

    mdict = mex3.maps[map_id]
    for k in ["T", "Wunpol"]:
        assert k in mdict
        assert mdict[k].compatible(m)
