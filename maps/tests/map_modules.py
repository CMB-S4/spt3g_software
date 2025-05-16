from spt3g import core, maps
import numpy as np

maplist = [
    maps.FlatSkyMap(300, 300, core.G3Units.arcmin, proj=maps.MapProjection.ProjZEA),
    maps.HealpixSkyMap(64),
]

for m in maplist:

    pipe = core.G3Pipeline()

    pipe.Add(core.G3InfiniteSource, type=core.G3FrameType.Observation, n=1)
    pipe.Add(maps.InjectMapStub, map_id="test_map", map_stub=m)

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

    pipe.Add(maps.ValidateMaps)
    pipe.Add(maps.SetPolConv, pol_conv=maps.MapPolConv.COSMO)
    pipe.Add(maps.ValidateMaps)

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

    coadder = maps.CoaddMaps(keep_outputs=True)
    pipe.Add(coadder)

    mex2 = maps.ExtractMaps(copy=True)
    pipe.Add(mex2)

    pipe.Add(maps.MakeMapsUnpolarized)

    mex3 = maps.ExtractMaps()
    pipe.Add(mex3)

    tmap = m.array_clone(np.random.randn(*m.shape))
    tmap.pol_type = maps.MapPolType.T
    tmap.weighted = False
    pipe.Add(maps.InjectMaps, map_id="test_map", maps_in=[tmap])

    mex4 = maps.ExtractMaps()
    pipe.Add(mex4)

    pipe.Run()

    # check that mex0 and mex1 are nearly (but not exactly) identical
    assert np.max(mex0.maps["test_map"]["T"] - mex1.maps["test_map"]["T"]) != 0
    assert np.allclose(mex0.maps["test_map"]["T"], mex1.maps["test_map"]["T"])

    # check that replication dropped the test map
    assert "test_map" not in mex2.maps
    assert "test_map" not in mex3.maps

    # check that injection added a new test map
    assert "test_map" in mex4.maps
    assert (mex4.maps["test_map"]["T"] == tmap).all()

    # check that coadd is correct
    assert "Coadd" in mex2.maps
    assert list(coadder.coadd_frame["InputMapIds"]) == ["test_map1", "test_map2"]
    assert (
        mex2.maps["Coadd"]["T"] ==
        mex2.maps["test_map1"]["T"] + mex2.maps["test_map2"]["T"]
    ).all()
    assert "Coadd" in mex3.maps
    assert (
        mex3.maps["Coadd"]["T"] ==
        mex3.maps["test_map1"]["T"] + mex3.maps["test_map2"]["T"]
    ).all()

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
