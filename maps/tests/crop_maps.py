#!/usr/bin/env python

import numpy as np
from spt3g import core, maps
from spt3g.maps import FlatSkyMap, MapProjection, MapCoordReference

# Map dimensions and the region of ones we'll paint in
NX, NY = 100, 80
RES = core.G3Units.arcmin

X1, X2 = 20, 45  # inclusive pixel bounds of the ones region
Y1, Y2 = 10, 35


def make_map(**kwargs):
    """Return a sparse FlatSkyMap with a rectangular region of ones."""
    m = FlatSkyMap(
        NX,
        NY,
        RES,
        proj=MapProjection.ProjZEA,
        coord_ref=MapCoordReference.Equatorial,
        **kwargs,
    )
    p = m[Y1 : Y2 + 1, X1 : X2 + 1]
    p += 1.0
    m[Y1 : Y2 + 1, X1 : X2 + 1] = p
    return m


# ------------------------------------------------------------------
# bounding_box tests
# ------------------------------------------------------------------

m = make_map()

# No pad: should recover exactly the region we painted
bbox = list(m.bounding_box())
expected = [Y1, Y2, X1, X2]
assert bbox == expected, f"bounding_box() returned {bbox}, expected {expected}"

# Pad of 3 pixels in both dimensions
PAD_PIX = 3
bbox_pad = list(m.bounding_box(x_pad=PAD_PIX, y_pad=PAD_PIX))
expected = [Y1 - PAD_PIX, Y2 + PAD_PIX, X1 - PAD_PIX, X2 + PAD_PIX]
assert bbox_pad == expected, f"bounding_box(pad=3) returned {bbox_pad}"

# Pad that would exceed the map edges should clamp to 0 / NX-1 / NY-1
bbox_clamp = list(m.bounding_box(x_pad=NX, y_pad=NY))
expected = [0, NY - 1, 0, NX - 1]
assert (
    bbox_clamp == expected
), f"bounding_box with large pad should clamp to map edges, got {bbox_clamp}"

# Empty map should return empty bounds
empty = FlatSkyMap(NX, NY, RES, proj=MapProjection.ProjZEA)
assert list(empty.bounding_box()) == [], "empty map should return empty bounding box"

print("bounding_box tests passed")


# ------------------------------------------------------------------
# crop tests
# ------------------------------------------------------------------

m = make_map()

# No pad: cropped map should be exactly (Y2-Y1+1) x (X2-X1+1)
cropped = m.crop()

expected_height = Y2 - Y1 + 1
expected_width = X2 - X1 + 1
expected_shape = (expected_height, expected_width)
assert (
    cropped.shape == expected_shape
), f"cropped shape {cropped.shape} != {expected_shape}"

# All values in the cropped map should be 1
arr = np.asarray(cropped)
assert np.all(arr == 1.0), "all cropped pixels should be 1"

# Crop with pad
PAD_PIX = 3
PAD_ANG = PAD_PIX * RES
cropped_pad = m.crop(pad=PAD_ANG)

padded_height = expected_height + 2 * PAD_PIX
padded_width = expected_width + 2 * PAD_PIX
expected_shape = (padded_height, padded_width)
assert (
    cropped_pad.shape == expected_shape
), f"padded crop shape {cropped_pad.shape} != {expected_shape}"

# Crop with pre-supplied bounds should give same result as crop without
cropped_b = m.crop_to(cropped)
assert np.array_equal(
    np.asarray(cropped_b), np.asarray(cropped)
), "crop_to() should match crop()"

# Empty map crop should return None
empty_cropped = empty.crop()
assert empty_cropped is None, "crop of empty map should return None"

print("crop tests passed")


# ------------------------------------------------------------------
# CropMaps pipeline module tests
# ------------------------------------------------------------------


def make_map_frame():
    stub = FlatSkyMap(
        NX,
        NY,
        RES,
        proj=MapProjection.ProjZEA,
        weighted=True,
        coord_ref=MapCoordReference.Equatorial,
    )
    fr = core.G3Frame(core.G3FrameType.Map)
    fr["Id"] = "test"

    tmap = make_map(weighted=True, pol_type=maps.MapPolType.T)
    fr["T"] = tmap

    wmap = maps.G3SkyMapWeights()
    tt = make_map(pol_type=maps.MapPolType.TT)
    wmap.TT = tt
    fr["Wunpol"] = wmap

    return fr


# Basic crop: frame maps should be trimmed to the observed region
fr = make_map_frame()
maps.CropMaps(fr)

expected = (Y2 - Y1 + 1, X2 - X1 + 1)
assert fr["T"].shape == expected, f"CropMaps T shape {fr['T'].shape}"
assert (
    fr["Wunpol"].TT.shape == expected
), f"CropMaps Wunpol.TT shape {fr['Wunpol'].TT.shape}"

# T and weight maps must be compatible after cropping
assert fr["T"].compatible(
    fr["Wunpol"].TT
), "T and Wunpol.TT must remain compatible after CropMaps"

# All values in cropped T should still be 1
assert np.all(np.asarray(fr["T"]) == 1.0), "CropMaps should preserve pixel values"

# With pad
fr_pad = make_map_frame()
maps.CropMaps(fr_pad, pad=PAD_ANG)

expected = (Y2 - Y1 + 1 + 2 * PAD_PIX, X2 - X1 + 1 + 2 * PAD_PIX)
assert fr_pad["T"].shape == expected, f"CropMaps with pad T shape {fr_pad['T'].shape}"
assert fr_pad["T"].compatible(
    fr_pad["Wunpol"].TT
), "T and Wunpol.TT must remain compatible after padded CropMaps"

print("CropMaps tests passed")
print("All tests passed")
