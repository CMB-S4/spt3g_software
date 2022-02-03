from spt3g import core, maps
import numpy as np

__all__ = [
    "CompactMaps",
    "RemoveWeights",
    "ApplyWeights",
    "FlattenPol",
    "MakeMapsPolarized",
    "MakeMapsUnpolarized",
    "ValidateMaps",
    "ExtractMaps",
    "InjectMapStub",
    "InjectMaps",
    "ReplicateMaps",
    "CoaddMaps",
]


@core.indexmod
def CompactMaps(frame, zero_nans=False):
    """
    Compact all maps in a frame to their default sparse representation.
    Optionally remove NaN values as well.  Removing NaN values will reduce
    memory use, but will remove the distinction in unweighted (or
    weight-removed) maps between unobserved regions and regions with zero
    temperature.
    """
    for s in ["T", "Q", "U", "Wunpol", "Wpol"]:
        if s in frame:
            m = frame.pop(s)
            m.compact(zero_nans=zero_nans)
            frame[s] = m
    return frame


@core.indexmod
def RemoveWeights(frame, zero_nans=False):
    """
    Remove weights from input maps.  If zero_nans is `True`, empty pixels are
    skipped and pixels with zero weight are set to 0 instead of NaN.  Operation
    is performed in place to minimize memory use.
    """
    if "Wpol" not in frame and "Wunpol" not in frame:
        return

    if not frame["T"].weighted:
        return frame

    tmap = frame.pop("T")

    if "Wpol" in frame:
        wmap = frame["Wpol"]
        qmap = frame.pop("Q")
        umap = frame.pop("U")
        maps.remove_weights(tmap, qmap, umap, wmap, zero_nans=zero_nans)
    else:
        wmap = frame["Wunpol"]
        maps.remove_weights_t(tmap, wmap, zero_nans=zero_nans)

    frame["T"] = tmap
    if "Wpol" in frame:
        frame["Q"] = qmap
        frame["U"] = umap

    return frame


@core.indexmod
def ApplyWeights(frame):
    """
    Apply weights to the input maps.  The operation is performed in place to
    minimize memory use.
    """
    if "Wpol" not in frame and "Wunpol" not in frame:
        return

    if frame["T"].weighted:
        return frame

    tmap = frame.pop("T")

    if "Wpol" in frame:
        wmap = frame["Wpol"]
        qmap = frame.pop("Q")
        umap = frame.pop("U")
        maps.apply_weights(tmap, qmap, umap, wmap)
    else:
        wmap = frame["Wunpol"]
        maps.apply_weights_t(tmap, wmap)

    frame["T"] = tmap
    if "Wpol" in frame:
        frame["Q"] = qmap
        frame["U"] = umap

    return frame


@core.indexmod
def FlattenPol(frame, invert=False):
    """
    For maps defined on the sphere the direction of the polarization angle is
    is defined relative to the direction of North.  When making maps we follow
    this definition.

    For any flat sky estimators, the polarization angle is defined relative to
    the vertical axis.  For some map projections the direction of north is not
    the same as the vertical axis.  This function applies a rotation to the Q
    and U values to switch the curved sky Q/U definition to the flat sky Q/U
    definition.

    If for whatever reason you want to reverse the process set the invert
    argument to True.
    """

    if "Q" not in frame or "U" not in frame:
        return

    if any(not isinstance(frame[k], maps.FlatSkyMap) for k in "QU"):
        return

    qmap, umap = frame.pop("Q"), frame.pop("U")

    if "Wpol" in frame:
        wmap = frame.pop("Wpol")
        maps.flatten_pol(qmap, umap, wmap, invert=invert)
        frame["Wpol"] = wmap
    else:
        maps.flatten_pol(qmap, umap, invert=invert)

    frame["Q"] = qmap
    frame["U"] = umap

    return frame


@core.indexmod
def MakeMapsPolarized(frame, pol_conv=maps.MapPolConv.IAU):
    """
    Converts individual unpolarized maps to polarized versions of the same map,
    with the given polarization convention

    This module is only a shim that creates null Q and U maps and populates
    a properly invertible Wpol array from the TT Wunpol weights.
    """
    if "Wunpol" not in frame:
        return

    wgt = frame["Wunpol"].TT
    del frame["Wunpol"]

    qmap = frame["T"].clone(False)
    qmap.pol_type = maps.MapPolType.Q
    frame["Q"] = qmap
    umap = frame["T"].clone(False)
    umap.pol_type = maps.MapPolType.U
    umap.pol_conv = pol_conv
    frame["U"] = umap
    mask = maps.get_mask_map(wgt)

    wgt_out = maps.G3SkyMapWeights(frame["T"], polarized=True)
    wgt_out.TT = wgt
    wgt_out.TQ = wgt.clone(False)
    wgt_out.TU = wgt.clone(False)
    wgt_out.QQ = mask
    wgt_out.QU = wgt.clone(False)
    wgt_out.UU = mask.clone(True)

    frame["Wpol"] = wgt_out

    return frame


@core.indexmod
def MakeMapsUnpolarized(frame):
    """
    Converts individual polarized maps to temperature-only versions of the same map.
    """
    if "Wpol" not in frame:
        return

    wgt = frame["Wpol"].TT
    del frame["Wpol"]
    del frame["Q"]
    del frame["U"]

    wgt_out = maps.G3SkyMapWeights(frame["T"], polarized=False)
    wgt_out.TT = wgt

    frame["Wunpol"] = wgt_out

    return frame


@core.indexmod
def ValidateMaps(frame, ignore_missing_weights=False):
    """
    Validate that the input map frame has all the necessary keys.

    If ignore_missing_weights is False (default), a warning is issued when the
    frame contains weighted Stokes maps without a weights map.  Set this option
    to True when feeding single bolometer map frames with common weights through
    a pipeline.
    """

    if isinstance(frame, core.G3Frame) and frame.type != core.G3FrameType.Map:
        return

    map_id = frame.get("Id", None)

    if "T" not in frame:
        core.log_fatal("Map frame %s: Missing T map" % map_id, unit="ValidateMaps")
    if ("Q" in frame and not "U" in frame) or ("U" in frame and not "Q" in frame):
        core.log_fatal("Map frame %s: Missing Q or U map" % map_id, unit="ValidateMaps")
    if "Wpol" in frame and "Wunpol" in frame:
        core.log_fatal(
            "Map frame %s: Found both polarized and unpolarized weights" % map_id,
            unit="ValidateMaps",
        )

    stub = frame["T"].clone(False)
    for k in ["T", "Q", "U", "Wpol", "Wunpol"]:
        if k not in frame:
            continue
        if not frame[k].compatible(stub):
            core.log_fatal(
                "Map frame %s: Map %s not compatible with T map" % (map_id, k),
                unit="ValidateMaps",
            )
        if k in "TQU":
            if k == "U" and frame[k].pol_conv is maps.MapPolConv.none:
                core.log_warn(
                    "Map frame %s: U map polarization convention not set" % map_id,
                    unit="ValidateMaps",
                )
            if frame[k].weighted and not ignore_missing_weights:
                if "Wpol" not in frame and "Wunpol" not in frame:
                    core.log_warn(
                        "Map frame %s: Missing weights" % map_id, unit="ValidateMaps"
                    )
                if k == "T" and "Q" not in frame and "Wunpol" not in frame:
                    core.log_warn(
                        "Map frame %s: Missing unpolarized weights" % map_id,
                        unit="ValidateMaps",
                    )
                if k in "QU" and "Wpol" not in frame:
                    core.log_warn(
                        "Map frame %s: Missing polarized weights" % map_id,
                        unit="ValidateMaps",
                    )
        else:
            if frame[k].polarized and ("Q" not in frame or "U" not in frame):
                core.log_fatal(
                    "Map frame %s: Found unpolarized maps with polarized weights"
                    % map_id,
                    unit="ValidateMaps",
                )
            elif not frame[k].polarized and ("Q" in frame or "U" in frame):
                core.log_fatal(
                    "Map frame %s: Found polarized maps with unpolarized weights"
                    % map_id,
                    unit="ValidateMaps",
                )


@core.indexmod
class ExtractMaps(object):
    """
    Cache maps that come through the pipeline. Initialize an instance of this
    module before adding to a pipeline..  Any maps that pass through the pipe
    are stored in the .maps attribute of the object after the pipeline is run.

    Arguments
    ---------
    map_id : string
        If supplied, select only map frames that match this ID.
    copy : bool
        If True, make a copy of the map on extraction.
    ignore_missing_weights : bool
        If False (default), a warning is issued when the frame contains weighted
        Stokes maps without a weights map.  Set this option to True when feeding
        single bolometer map frames with common weights through a pipeline.
    """

    def __init__(self, map_id=None, copy=False, ignore_missing_weights=False):
        self.map_id = map_id
        self.copy_ = copy
        self.ignore_missing_weights = ignore_missing_weights
        self.maps = {}

    def __call__(self, frame):
        if frame.type != core.G3FrameType.Map:
            return
        if self.map_id and frame["Id"] != self.map_id:
            return

        ValidateMaps(frame, ignore_missing_weights=self.ignore_missing_weights)

        mid = frame["Id"]
        mdict = {}
        for k in ["T", "Q", "U", "Wpol", "Wunpol"]:
            if k not in frame:
                continue
            mdict[k] = frame[k] if not self.copy_ else frame[k].copy()

        if mid not in self.maps:
            self.maps[mid] = mdict
            return

        if isinstance(self.maps[mid], dict):
            self.maps[mid] = [self.maps[mid], mdict]
            return

        self.maps[mid].append(mdict)


@core.indexmod
class InjectMapStub(object):
    """
    Inject a new map frame from a map stub.

    Arguments
    ---------
    map_id : string
        Id to assign to the new map frame
    map_stub : G3SkyMap instance
        Map stub from which to clone the Stokes maps and weights.
    polarized : bool
        If True, add Q and U maps to stub frame, and ensure that weights are
        polarized.  Otherwise, only a T map is created.
    weighted : bool
        If True, add weights to the stub frame.
    pol_conv : MapPolConv instance
        Polarization convention to use.
    """

    def __init__(
        self,
        map_id,
        map_stub,
        polarized=True,
        weighted=True,
        pol_conv=maps.MapPolConv.IAU,
    ):
        self.map_frame = core.G3Frame(core.G3FrameType.Map)
        self.map_frame["Id"] = map_id

        map_stub = map_stub.clone(False)
        map_stub.weighted = weighted
        map_stub.pol_conv = pol_conv

        T = map_stub.clone(False)
        T.pol_type = maps.MapPolType.T
        self.map_frame["T"] = T
        if polarized:
            Q = map_stub.clone(False)
            Q.pol_type = maps.MapPolType.Q
            self.map_frame["Q"] = Q
            U = map_stub.clone(False)
            U.pol_type = maps.MapPolType.U
            self.map_frame["U"] = U
        if weighted:
            W = maps.G3SkyMapWeights(map_stub, polarized)
            self.map_frame["Wpol" if polarized else "Wunpol"] = W

    def __call__(self, frame):
        if self.map_frame is None:
            return

        map_frame = self.map_frame
        self.map_frame = None
        return [map_frame, frame]


@core.indexmod
class InjectMaps(object):
    """
    Inject a set of maps into a new map frame.

    Arguments
    ---------
    map_id : string
        Id to assign to the new map frame
    maps_in : list or dict
        Maps to add to the frame.  If a list, contains Stokes maps with valid
        pol_type and weights.  If a dict, contains Stokes and weights maps keyed
        by the standard map frame names.
    ignore_missing_weights [False] : bool
        Skip warning about missing weights.  Useful for masks.
    """

    def __init__(self, map_id, maps_in, ignore_missing_weights=False):
        self.map_frame = core.G3Frame(core.G3FrameType.Map)
        self.map_frame["Id"] = map_id

        if isinstance(maps_in, list):
            for m in maps_in:
                if isinstance(m, maps.G3SkyMap):
                    k = str(m.pol_type)
                    if k not in "TQU":
                        raise ValueError("Input map has invalid pol_type %s" % k)
                    self.map_frame[k] = m
                elif isinstance(m, maps.G3SkyMapWeights):
                    self.map_frame["Wpol" if m.polarized else "Wunpol"] = m
                else:
                    raise TypeError("maps_in must be G3SkyMap or G3SkyMapWeights")

        elif isinstance(maps_in, dict):
            for k, m in maps_in.items():
                if k not in ["T", "Q", "U", "Wpol", "Wunpol"]:
                    continue
                self.map_frame[k] = m

        else:
            raise TypeError("maps_in must be a list or dict")

        ValidateMaps(self.map_frame, ignore_missing_weights=ignore_missing_weights)

    def __call__(self, frame):
        if self.map_frame is None:
            return

        map_frame = self.map_frame
        self.map_frame = None
        return [map_frame, frame]


@core.indexmod
def ReplicateMaps(frame, input_map_id, output_map_ids, copy_weights=False):
    """
    Clone the input map frame with Id input_map_id into new stub frames, one for
    each Id listed in output_map_ids.

    Arguments
    ---------
    input_map_id : string
        ID of the map frame to replicate.  The input frame is discarded after
        replication.
    output_map_ids : list of strings
        List of IDs to assign to replicated map frames.
    copy_weights : bool
        If False, only the first output frame in the list includes a weights key
        (Wpol or Wunpol).  If True, all output frames include a weights key.
    """

    if frame.type != core.G3FrameType.Map:
        return

    if frame["Id"] != input_map_id:
        return

    ValidateMaps(frame)

    frames = []

    first = True
    for oid in output_map_ids:
        fr = core.G3Frame(core.G3FrameType.Map)
        fr["Id"] = oid
        if copy_weights or first:
            map_keys = ["T", "Q", "U", "Wpol", "Wunpol"]
            first = False
        else:
            map_keys = ["T", "Q", "U"]

        for k in map_keys:
            if k not in frame:
                continue
            fr[k] = frame[k].clone(False)

        frames.append(fr)

    return frames


@core.indexmod
class CoaddMaps(object):
    """
    Coadd maps and weights.

    Arguments
    ---------
    map_ids : list of str
        List of map Id's to include in the coadd.  If None, any maps
        in the pipeline are included.
    output_map_id : str
        Id to assign to the output frame.
    ignore_missing_weights : bool
        If False (default), a warning is issued when the frame contains weighted
        Stokes maps without a weights map.  Set this option to True when feeding
        single bolometer map frames with common weights through a pipeline.
    """

    def __init__(self, map_ids=None, output_map_id=None, ignore_missing_weights=False):
        self.coadd_frame = core.G3Frame(core.G3FrameType.Map)
        self.coadd_frame["Id"] = output_map_id
        if isinstance(map_ids, str):
            map_ids = [map_ids]
        self.map_ids = map_ids
        self.ignore_missing_weights = ignore_missing_weights

    def __call__(self, frame):

        if frame.type == core.G3FrameType.EndProcessing:
            coadd = self.coadd_frame
            self.coadd_frame = None
            return [coadd, frame]

        if "Id" not in frame:
            return

        if self.map_ids is not None and frame["Id"] not in self.map_ids:
            return

        ValidateMaps(frame, ignore_missing_weights=self.ignore_missing_weights)
        input_weighted = True
        if not frame["T"].weighted:
            input_weighted = False
            ApplyWeights(frame)

        for key in ["T", "Q", "U", "Wpol", "Wunpol"]:
            if key not in frame:
                continue
            if key not in self.coadd_frame:
                self.coadd_frame[key] = frame[key].clone(False)
            m = self.coadd_frame.pop(key)
            m += frame[key]
            self.coadd_frame[key] = m

        if not input_weighted:
            RemoveWeights(frame)
