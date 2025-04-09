#!/usr/bin/env python

from spt3g import core
import numpy as np

dtypes = [np.int32, np.int64, np.float32, np.float64]
start = core.G3Time.Now()
stop = start + 5 * core.G3Units.s
ndet = 20
nsamp = 1200
keys = [str(x).zfill(3) for x in range(ndet)]

for bit_depth in [24, 32]:
    for dtype in dtypes:
        print("bit depth", bit_depth, "dtype", dtype)
        dat = np.random.normal(size=(ndet, nsamp), scale=2**28, loc=0).astype(dtype)
        if dtype in [np.float32, np.float64]:
            dat[:, 5] = np.nan
        try:
            tsm = core.G3TimestreamMap(
                keys, dat, start=start, stop=stop, compression_level=5, bit_depth=bit_depth
            )
        except Exception as e:
            if "32-bit compression is not supported" in str(e):
                continue
            raise

        with core.G3Writer("test.g3") as w:
            fr = core.G3Frame()
            fr["tsm"] = tsm
            w(fr)

        tsm2 = list(core.G3File("test.g3"))[0]["tsm"]

        assert tsm.start == tsm2.start
        assert tsm.stop == tsm2.stop
        assert tsm.compression_level == tsm2.compression_level
        assert tsm.bit_depth == tsm2.bit_depth

        if bit_depth == 24 and tsm.dtype == np.float64:
            assert tsm2.dtype == np.float32
        elif bit_depth == 24 and tsm.dtype == np.int64:
            assert tsm2.dtype == np.int32
        else:
            assert tsm.dtype == tsm2.dtype

        if dtype in [np.float32, np.float64]:
            assert np.isnan(np.asarray(tsm2)[:, 5]).all()

        def mask32(v):
            return (np.asarray(v)[:, 6:]).astype(np.int32)

        tsm_trunc = mask32(tsm)
        if bit_depth == 24:
            tsm_trunc = ((tsm_trunc & 0x00FFFFFF) << 8) >> 8
        np.testing.assert_array_equal(tsm_trunc, mask32(tsm2))
