from spt3g import core, maps
import numpy as np

__all__ = [
    "CompactMaps",
    "RemoveWeights",
    "ApplyWeights",
    "SetPolConv",
    "FlattenPol",
    "MakeMapsPolarized",
    "MakeMapsUnpolarized",
    "ValidateMaps",
    "ExtractMaps",
    "InjectMapStub",
    "InjectMaps",
    "ReplicateMaps",
    "CoaddMaps",
    "ReprojectMaps",
    "coadd_map_files",
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
    for s in ["T", "Q", "U", "Wunpol", "Wpol", "H"]:
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
    ValidateMaps(frame)

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
    ValidateMaps(frame)

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
def SetPolConv(frame, pol_conv=maps.MapPolConv.IAU):
    """
    Set or change the polarization convention of the input polarized map frame.
    If switching between IAU and COSMO conventions, flip the sign of the U map
    and the TU and QU weights.  Otherwise, just set the polarization convention
    for all maps and weights.
    """
    if not ("Q" in frame and "U" in frame):
        # only polarized frames
        return frame

    if pol_conv == maps.MapPolConv.none or pol_conv is None:
        raise ValueError("Polarized maps must have pol_conv set to IAU or COSMO")

    tmap = frame.pop("T")
    tmap.pol_conv = pol_conv

    qmap = frame.pop("Q")
    qmap.pol_conv = pol_conv

    umap = frame.pop("U")
    flip = umap.polarized and umap.pol_conv != pol_conv
    umap.pol_conv = pol_conv

    wmap = None
    if "Wpol" in frame:
        wmap = frame.pop("Wpol")
        for k in wmap.keys():
            wmap[k].pol_conv = pol_conv

    # flip sign if switching conventions
    if flip:
        umap *= -1.0
        if wmap is not None:
            wmap.TU *= -1.0
            wmap.QU *= -1.0

    frame["T"] = tmap
    frame["Q"] = qmap
    frame["U"] = umap
    if wmap is not None:
        frame["Wpol"] = wmap

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

    ValidateMaps(frame)
    tmap, qmap, umap = frame.pop("T"), frame.pop("Q"), frame.pop("U")

    if "Wpol" in frame:
        wmap = frame.pop("Wpol")
        maps.flatten_pol(qmap, umap, wmap, invert=invert)
        frame["Wpol"] = wmap
    else:
        maps.flatten_pol(qmap, umap, invert=invert)

    tmap.flat_pol = not invert
    
    frame["T"] = tmap
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
    wgt.pol_conv = pol_conv
    del frame["Wunpol"]

    qmap = frame["T"].clone(False)
    qmap.pol_type = maps.MapPolType.Q
    qmap.pol_conv = pol_conv
    frame["Q"] = qmap
    umap = frame["T"].clone(False)
    umap.pol_type = maps.MapPolType.U
    umap.pol_conv = pol_conv
    frame["U"] = umap
    mask = wgt.to_mask().to_map()

    wgt_out = maps.G3SkyMapWeights()
    wgt_out.TT = wgt
    wgt_out.TQ = wgt.clone(False)
    wgt_out.TU = wgt.clone(False)
    wgt_out.QQ = mask
    wgt_out.QU = wgt.clone(False)
    wgt_out.UU = mask.clone(True)

    for k in wgt_out.keys():
        wgt_out[k].pol_type = getattr(maps.MapPolType, k)

    frame["Wpol"] = wgt_out

    return frame


@core.indexmod
def MakeMapsUnpolarized(frame):
    """
    Converts individual polarized maps to temperature-only versions of the same map.
    """
    if "Wpol" not in frame:
        return

    tmap = frame.pop("T")
    tmap.pol_conv = maps.MapPolConv.none
    frame["T"] = tmap

    wgt = frame.pop("Wpol").TT
    wgt.pol_conv = maps.MapPolConv.none
    wgt_out = maps.G3SkyMapWeights()
    wgt_out.TT = wgt
    frame["Wunpol"] = wgt_out

    del frame["Q"]
    del frame["U"]

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
        if "H" in frame:
            return
        core.log_fatal("Map frame %s: Missing T map" % map_id, unit="ValidateMaps")
    if ("Q" in frame and not "U" in frame) or ("U" in frame and not "Q" in frame):
        core.log_fatal("Map frame %s: Missing Q or U map" % map_id, unit="ValidateMaps")
    if "Wpol" in frame and "Wunpol" in frame:
        core.log_fatal(
            "Map frame %s: Found both polarized and unpolarized weights" % map_id,
            unit="ValidateMaps",
        )

    check_weights = False

    stub = frame["T"].clone(False)
    for k in ["T", "Q", "U", "Wpol", "Wunpol", "H"]:
        if k not in frame:
            continue
        if not frame[k].compatible(stub):
            core.log_fatal(
                "Map frame %s: Map %s not compatible with T map" % (map_id, k),
                unit="ValidateMaps",
            )

        if k in ["Wpol", "Wunpol"]:
            if frame[k].TT.pol_type == maps.MapPolType.TT:
                continue
            # set weights polarization properties
            w = frame.pop(k)
            for wk in w.keys():
                w[wk].pol_type = getattr(maps.MapPolType, wk)
                if k == "Wpol":
                    w[wk].pol_conv = frame["U"].pol_conv
            frame[k] = w

        if k in "TQU":
            if k == "U":
                if isinstance(frame[k], maps.FlatSkyMap) and (
                    frame[k].flat_pol != frame["Q"].flat_pol
                ):
                    core.log_fatal(
                        "Map frame %s: Q and U maps have different flat_pol" % map_id,
                        unit="ValidateMaps",
                    )
            if k in "QU":
                if not frame[k].polarized:
                    core.log_warn(
                        "Map frame %s: %s map polarization convention not set" % (map_id, k),
                        unit="ValidateMaps",
                    )
            if frame[k].weighted and not ignore_missing_weights:
                if "Wpol" not in frame and "Wunpol" not in frame:
                    if not check_weights:
                        core.log_warn(
                            "Map frame %s: Missing weights" % map_id, unit="ValidateMaps"
                        )
                        check_weights = True
                else:
                    if k == "T" and "Q" not in frame and "Wunpol" not in frame:
                        core.log_warn(
                            "Map frame %s: Missing unpolarized weights" % map_id,
                            unit="ValidateMaps",
                        )
                    elif k == "Q" and "Wpol" not in frame:
                        core.log_warn(
                            "Map frame %s: Missing polarized weights" % map_id,
                            unit="ValidateMaps",
                        )

        elif k == "H":
            continue

        elif frame[k].polarized:
            if "Q" not in frame or "U" not in frame:
                core.log_fatal(
                    "Map frame %s: Found unpolarized maps with polarized weights" % map_id,
                    unit="ValidateMaps",
                )
            if frame[k].pol_conv != frame["U"].pol_conv:
                core.log_fatal(
                    "Map frame %s: %s and U maps have different pol_conv" % (map_id, k),
                    unit="ValidateMaps",
                )
            if isinstance(frame[k].QQ, maps.FlatSkyMap):
                if frame[k].flat_pol != frame["Q"].flat_pol:
                    core.log_fatal(
                        "Map frame %s: %s and U maps have different flat_pol" % (map_id, k),
                        unit="ValidateMaps",
                    )

        elif "Q" in frame or "U" in frame:
            core.log_fatal(
                "Map frame %s: Found polarized maps with unpolarized weights" % map_id,
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
        for k in ["T", "Q", "U", "Wpol", "Wunpol", "H"]:
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
        Map stub from which to clone the Stokes maps and weights.  If the
        `weighted` attribute of the stub is True, the output frame will include
        weights.  If the `pol_conv` attribute of the stub is not None, the
        output frame will include Q and U maps (and polarized weights).
    """

    def __init__(self, map_id, map_stub):
        self.map_frame = core.G3Frame(core.G3FrameType.Map)
        self.map_frame["Id"] = map_id

        map_stub = map_stub.clone(False)
        weighted = map_stub.weighted
        polarized = map_stub.polarized

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
            W = maps.G3SkyMapWeights(map_stub)
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
                if k not in ["T", "Q", "U", "Wpol", "Wunpol", "H"]:
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
            map_keys = ["T", "Q", "U", "Wpol", "Wunpol", "H"]
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
    Coadd maps and weights, optionally collating by map Id.  This class can be
    used as an argument to ``pipe.Add()`` as a standard pipeline module, or
    instantiated as a standalone instance.  In the latter case, the object is
    treated as a callable for input map frames, and the ``coadd_frame`` or
    ``coadd_frames`` attribute contains the running coadd(s).

    The output coadd frames contain two additional keys: ``'InputMapIds'`` and
    ``'InputFiles'``, which are both lists of unique map Id's and filenames that
    are associated with the frames that contribute to each coadd.  When one
    coadd is added to another, these keys are updated recursively, so that the
    resulting coadd includes the Id's and filenames that contributed to both
    itself and the other coadd.  The list of filenames can be populated by
    combining this module with a G3Reader whose ``track_filename`` option is set
    to True; however, this feature is fragile and may not work as expected with
    complex pipelines.

    Attributes
    ----------
    coadd_frame : G3Frame
        Output coadd map frame, also injected into the pipeline on
        EndProcessing.  This attribute is only populated if the ``collate``
        option is set to False.
    coadd_frames : dict of G3Frames
        Output coadd map frames, keyed by input map Id.  Each frame is also
        injected into the pipeline on EndProcessing.  This attribute is only
        populated if the ``collate`` option is set to True.

    Methods
    -------
    get_map_id :
        Takes a map frame as an argument and returns an identifier string for
        the coadd to which it should be added, or None if the map should be
        ignored.  This method can be modified by subclassing the CoaddMaps
        module.

    Arguments
    ---------
    map_ids : list of str
        List of map Id's to include in the coadd(s).  If None, any maps in the
        pipeline are included.  Otherwise, the output of the ``get_map_ids``
        method is compared with this list, and the input frame is discarded if
        no match is found.
    output_map_id : str
        Id to assign to the output frame.  If ``collate`` is True, this argument
        is required and treated as a prefix to which each input map Id is
        appended.
    collate : bool
        If True, coadd unique map Id's into separate output map frames.
    weighted : bool
        If True (default), ensure that maps have had weights applied before
        coadding.  Otherwise, coadd maps without checking the weights.
    ignore_missing_weights : bool
        If False (default), a warning is issued when the frame contains weighted
        Stokes maps without a weights map.  Set this option to True when feeding
        single bolometer map frames with common weights through a pipeline.
    drop_input_frames : bool
        If True, drop input map frames from the pipeline that are included in
        any coadds.
    record_obs_id : bool
        If True, include source name and observation ID info in the output coadd
        frame ``InputMapIds`` key, along with the map ID for each input frame.
        If False, only the map frame ID is included.
    """

    def __init__(
        self,
        map_ids=None,
        output_map_id="Coadd",
        collate=False,
        weighted=True,
        ignore_missing_weights=False,
        drop_input_frames=False,
        record_obs_id=False,
    ):
        if isinstance(map_ids, str):
            map_ids = [map_ids]
        self.map_ids = map_ids
        self.collate = collate
        if self.collate:
            self.coadd_frames = dict()
            self.output_map_id = output_map_id
        else:
            self.coadd_frame = core.G3Frame(core.G3FrameType.Map)
            self.coadd_frame["Id"] = output_map_id
        self.weighted = weighted
        self.ignore_missing_weights = ignore_missing_weights
        self.drop_input_frames = drop_input_frames
        self.obs_id = None if record_obs_id else False

    def get_map_id(self, frame):
        """
        Return Id associated with the input frame

        By default, this returns the "Id" entry in the frame, None if not found.
        Subclass the CoaddMaps structure to override the default behavior.

        Arguments
        ---------
        frame : core.G3Frame instance of type G3FrameType.Map
            Candidate map frame to include in a coadd.

        Returns
        -------
        map_id : str, None or False
            A string if the frame is to be included in a coadd, None if the
            frame is to be passed on to downstream pipeline modules, or False if
            the frame is to be dropped from the pipeline altogether.
        """
        return frame.get("Id", None)

    def __call__(self, frame):

        if isinstance(frame, core.G3Frame) and frame.type == core.G3FrameType.EndProcessing:
            if self.collate:
                return list(self.coadd_frames.values()) + [frame]
            return [self.coadd_frame, frame]

        if self.obs_id is not False and "SourceName" in frame:
            self.obs_id = "{}/{}".format(
                frame["SourceName"], frame["ObservationID"]
            )

        if isinstance(frame, core.G3Frame) and frame.type != core.G3FrameType.Map:
            return

        map_id = self.get_map_id(frame)
        if map_id is None or map_id is False:
            return map_id

        if self.map_ids is not None and map_id not in self.map_ids:
            return

        ValidateMaps(frame, ignore_missing_weights=self.ignore_missing_weights)
        if self.weighted:
            ApplyWeights(frame)

        if self.collate:
            if map_id not in self.coadd_frames:
                fr = core.G3Frame(core.G3FrameType.Map)
                fr["Id"] = self.output_map_id + map_id
                self.coadd_frames[map_id] = fr
            cfr = self.coadd_frames[map_id]
        else:
            cfr = self.coadd_frame

        if "InputMapIds" in cfr:
            map_ids = list(cfr.pop("InputMapIds"))
        else:
            map_ids = []
        if "InputMapIds" in frame:
            # allow for recursive coadds
            map_ids += [i for i in frame["InputMapIds"] if i not in map_ids]
        elif frame.get("Id", None):
            mid = frame["Id"]
            if self.obs_id:
                mid = "{}/{}".format(self.obs_id, mid)
            if mid not in map_ids:
                map_ids += [mid]
        if len(map_ids):
            cfr["InputMapIds"] = core.G3VectorString(map_ids)

        if "InputFiles" in cfr:
            input_files = list(cfr.pop("InputFiles"))
        else:
            input_files = []
        if "InputFiles" in frame:
            # allow for recursive coadds
            input_files += [
                f for f in frame["InputFiles"] if f not in input_files
            ]
        elif getattr(frame, "_filename", None):
            if frame._filename not in input_files:
                input_files += [frame._filename]
        if len(input_files):
            cfr["InputFiles"] = core.G3VectorString(input_files)

        for key in ["T", "Q", "U", "Wpol", "Wunpol", "H"]:
            if key not in frame:
                continue
            if key not in cfr:
                cfr[key] = frame[key].clone(False)
            m = cfr.pop(key)
            m += frame[key]
            cfr[key] = m

        if self.drop_input_frames:
            return False


@core.usefulfunc
def coadd_map_files(
    input_files,
    output_file=None,
    coadder=None,
    map_ids=None,
    output_map_id="Coadd",
    collate=False,
    weighted=True,
    record_obs_id=False,
):
    """
    Coadd map files, optionally collating map Id's into separate frames.

    Arguments
    ---------
    input_files : list of str
        List of input files to feed through the pipeline.
    output_file : str
        Output G3 filename.  If not supplied, the output frames are
        returned without saving to disk.
    coadder : CoaddMaps instance
        If set, use this instantiated module in the coadding pipeline.
        In this case, all other keyword arguments below are ignored.
    map_ids : list of str
        A list of map Id's to include in the coadd(s).
    output_map_id : str
        Id to use for the output map frame.  If ``collate`` is True,
        this is the prefix applied to each output frame, with the
        input map Id appended to it.
    collate : bool
        If True, coadd individual map Id's into separate map frames.
        Otherwise, all map frames are coadded into one output frame.
    weighted : bool
        If True, ensure that weights have been applied before coadding.
        Otherwise, the input maps are coadded as they are.
    record_obs_id : bool
        If True, include source name and observation ID info in the output coadd
        frame ``InputMapIds`` key, along with the map ID for each input frame.
        If False, only the map frame ID is included.

    Returns
    -------
    maps : G3Frame or dict of G3Frames
        If ``collate`` is True, returns a dictionary of map frames
        keyed by Id.  Otherwise, returns a single map frame.
    """

    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=input_files, track_filename=True)

    # build coadds
    if coadder is None:
        coadder = CoaddMaps(
            map_ids=map_ids,
            output_map_id=output_map_id,
            collate=collate,
            weighted=weighted,
            drop_input_frames=True,
            record_obs_id=record_obs_id,
        )
    pipe.Add(coadder)

    # drop metadata frames
    pipe.Add(lambda fr: fr.type == core.G3FrameType.Map)

    if output_file:
        pipe.Add(core.G3Writer, filename=output_file)
    pipe.Run()

    if hasattr(coadder, 'coadd_frames'):
        return coadder.coadd_frames
    return coadder.coadd_frame


@core.indexmod
class ReprojectMaps(object):
    """
    Reproject a map frame into a different projection.  Original data are
    dropped and replaced by reprojected maps in the input frames.  Maps can be
    changed between flat sky and healpix pixelizations, rotated between
    Equatorial and Galactic coordinates, and/or change polarization convention
    between COSMO and IAU, by setting the appropriate attributes of the input
    and stub maps.  Attributes not defined in the stub map are assumed to be
    that of the input map.  NB: coordinate rotation of polarized maps is not
    currently implemented.

    Arguments
    ---------
    map_stub : G3SkyMap object
        A stub (empty) sky map object to be used to construct the output maps.
        Can be a HealpixSkyMap or FlatSkyMap object.  Setting the ``pol_conv``
        and/or ``coord_ref`` attributes to values that differ from those of the
        input maps will result in output maps whose polarization convention
        and/or reference coordinate system have been changed.
    rebin : int
        If supplied and >1, subdivide the output pixel by n x n with each
        sub-pixel taking on the input map values at pixel center (with interp or
        nearest neighbor). The output pixel takes on the average of the
        sub-pixel values.  In the case that the input map has higher resolution
        than the output map (and that the input map is not low-pass filtered to
        remove information above the Nyquist freq. of the output map pixel),
        this reduces aliasing compared with direct sampling. But there would
        still be aliased power from the input map from freq above the ouput map
        pixel's Nyquist.
    interp : bool
        If True, use bilinear interpolation to extract values from the input
        map.  Otherwise, the nearest-neighbor value is used.
    """

    def __init__(self, map_stub=None, rebin=1, interp=False):
        assert map_stub is not None, "map_stub argument required"
        self.stub = map_stub
        self.rebin = rebin
        self.interp = interp

    def __call__(self, frame):

        if isinstance(frame, core.G3Frame) and frame.type != core.G3FrameType.Map:
            return

        if "Q" in frame and self.stub.coord_ref != frame["Q"].coord_ref:
            raise RuntimeError(
                "Coordinate rotation of polarized maps is not implemented"
            )

        if "U" in frame and not self.stub.polarized:
            self.stub.pol_conv = frame["U"].pol_conv

        for key in ["T", "Q", "U", "Wpol", "Wunpol", "H"]:

            if key not in frame:
                continue

            m = frame.pop(key)

            if key in "TQUH":
                mnew = self.stub.clone(False)
                maps.reproj_map(m, mnew, rebin=self.rebin, interp=self.interp)

            elif key in ["Wpol", "Wunpol"]:
                mnew = maps.G3SkyMapWeights(self.stub)
                for wkey in mnew.keys():
                    maps.reproj_map(
                        m[wkey], mnew[wkey], rebin=self.rebin, interp=self.interp
                    )

            frame[key] = mnew

        return frame
