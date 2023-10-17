from spt3g import core
from spt3g.core import G3TimestreamQuat
from spt3g.maps import MapCoordReference
from spt3g.maps import get_origin_rotator_timestream, get_boresight_rotator_timestream
from spt3g.maps import get_fk5_j2000_to_gal_quat


__all__ = [
    "FillCoordTransRotations",
    "EquatorialToGalacticTransRotations",
    "AddLocalTransRotations",
]


@core.indexmod
def FillCoordTransRotations(
    frame,
    transform_store_key="OnlineRaDecRotation",
    end_coord_sys=MapCoordReference.Equatorial,
    do_bad_transform=False,
    bs_az_key="RawBoresightAz",
    bs_el_key="RawBoresightEl",
    bs_ra_key="OnlineBoresightRa",
    bs_dec_key="OnlineBoresightDec",
    offset_az_key="OffsetBoresightAz",
    offset_el_key="OffsetBoresightEl",
    offset_ra_key="OnlineOffsetRa",
    offset_dec_key="OnlineOffsetDec",
):
    """
    Calculates the rotation quaternions that take the point (1,0,0) (so az=el=0)
    in local coordinates to the coordinates specified by end_coord_sys and
    stores them in transform_store_key.  This encodes the boresight pointing and
    any rotations about this boresight pointing due to coordinate system
    changes, az/el bearing tilt, etc.

    Arguments
    ---------
    transform_store_key : string
        The key where the output transformation quaternion will be stored.
        If already present in the frame, this calculation will be skipped.
    end_coord_sys : MapCoordReference
        If Local, the transformation is computed using the negative of the
        detector delta angle.  Otherwise the detector angle is not inverted.
    do_bad_transform : bool
        If end_coord_sys is not Local and this argument is True, the offset
        keys are ignored and the coordinate transformation does not take
        into account rotation about the boresight.
    bs_az_key, bs_el_key : string
        Boresight coordinates in the local coordinate system.  If end_coord_sys
        is Local, only these two keys are required.
    bs_ra_key, bs_dec_key : string
        Boresight coordinates in the output coordinate system.  If
        do_bad_transform is True, only these two keys and the previous two keys
        are required.
    offset_az_key, offset_el_key, offset_ra_key, offset_dec_key : string
        Local and output coordinates computed at a small offset from boresight.
        These keys are required if do_bad_transform is False, and will be used
        to account for any rotation about boresight in the coordinate
        transformation.
    """
    if frame.type != core.G3FrameType.Scan:
        return

    if transform_store_key in frame:
        core.log_debug("Transform already computed, skipping")
        return

    trans = G3TimestreamQuat()
    if end_coord_sys == MapCoordReference.Local:
        trans = get_origin_rotator_timestream(
            frame[bs_az_key], frame[bs_el_key], end_coord_sys
        )
    else:
        if do_bad_transform:
            core.log_debug("You are doing the old calculation for pointing")
            trans = get_origin_rotator_timestream(
                frame[bs_ra_key], frame[bs_dec_key], end_coord_sys
            )
        else:
            trans = get_boresight_rotator_timestream(
                frame[bs_az_key],
                frame[bs_el_key],
                frame[bs_ra_key],
                frame[bs_dec_key],
                frame[offset_az_key],
                frame[offset_el_key],
                frame[offset_ra_key],
                frame[offset_dec_key],
            )
    frame[transform_store_key] = trans


@core.indexmod
def EquatorialToGalacticTransRotations(
    frame, eq_trans_key="OnlineRaDecRotation", out_key="OnlineGalacticRotation"
):
    """
    Takes a quaternion vector specifying the rotation to FK5 (Equatorial)
    boresight and converts it into a rotation to Galactic J2000 boresight
    coordinates.

    Use this to convert the output of FillCoordTransRotations to Galactic
    coordinates.
    """

    if frame.type != core.G3FrameType.Scan:
        return
    gal_trans = get_fk5_j2000_to_gal_quat() * frame[eq_trans_key]
    frame[out_key] = gal_trans


@core.indexmod
def AddLocalTransRotations(
    frame, az_key="RawBoresightAz", el_key="RawBoresightEl", out_key="RawAzElRotation"
):
    """
    Creates the transform for boresight pointing for az el based maps.  This is
    equivalent to FillCoordTransRotations with end_coord_sys in Local
    coordinates.

    Right now it's using the coordinate system where delta = -el because of
    implementation details.  All of the maps generated with this will be upside
    down in el.
    """

    if frame.type != core.G3FrameType.Scan:
        return
    local_trans = get_origin_rotator_timestream(
        frame[az_key], frame[el_key], MapCoordReference.Local
    )
    frame[out_key] = local_trans

