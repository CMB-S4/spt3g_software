from spt3g import core

from spt3g.core import G3VectorQuat

from spt3g.coordinateutils.azel import convert_azel_to_radec

from spt3g.coordinateutils import MapCoordReference
from spt3g.coordinateutils import create_det_az_el_trans, create_lazy_det_ra_dec_trans
from spt3g.coordinateutils import create_det_ra_dec_trans, convert_ra_dec_trans_to_gal

from copy import copy
import numpy as np

@core.indexmod
def FillCoordTransRotations(
        frame, 
        transform_store_key,
        end_coord_sys = MapCoordReference.Equatorial,
        
        do_bad_transform = False,
        
        bs_az_key = 'RawBoresightAz', bs_el_key='RawBoresightEl',
        bs_ra_key = 'OnlineBoresightRa', bs_dec_key='OnlineBoresightDec',
        offset_az_key='OffsetBoresightAz', offset_el_key='OffsetBoresightEl', 
        offset_ra_key='OnlineOffsetRa', offset_dec_key='OnlineOffsetDec'):
    '''
    Calculates the rotation quaternions that take the point (1,0,0) (so az=el=0)
    in local coordinates to the coordinates specified by end_coord_sys 
    and stores them in transform_store_key.  This encodes the boresight pointing and
    any rotations about this boresight pointing due to coordinate system changes, az/el bearing
    tilt, etc.

    If do_bad_transform:
      It uses one set of points to calculate this rotations:
      (bs_az_key, bs_el_key) (bs_ra_key, bs_dec_key) 
    Else:
      It uses two sets of points to calculate this rotations:
      (bs_az_key, bs_el_key) (bs_ra_key, bs_dec_key) 
      (offset_az_key, offset_el_key) (offset_ra_key, offset_dec_key) 

    If you do a bad transform it will not properly calculate the rotation around the
    boresight axis.
    '''
    if frame.type != core.G3FrameType.Scan:
        return

    if transform_store_key in frame:
        core.log_debug("Transform already computed, skipping")
        return

    trans = G3VectorQuat()
    if end_coord_sys == MapCoordReference.Equatorial:
        if do_bad_transform:
            core.log_debug("You are doing the old calculation for pointing")
            create_lazy_det_ra_dec_trans(frame[bs_ra_key], frame[bs_dec_key], trans)
        else:
            create_det_ra_dec_trans(
                frame[bs_az_key], frame[bs_el_key],
                frame[bs_ra_key], frame[bs_dec_key],
                frame[offset_az_key], frame[offset_el_key],
                frame[offset_ra_key], frame[offset_dec_key],
                trans)
    elif end_coord_sys == MapCoordReference.Local:
        create_det_az_el_trans(frame[bs_az_key], frame[bs_el_key], trans)
    else:
        core.logfatal("To do the galactic convert the equatorial transform to galactic later "+
                      "using convert_ra_dec_trans_to_gal")
    frame[transform_store_key] = trans

def AddGalTrans(frame, eq_trans_key, out_key):
    '''
    Takes a quaternion vector specifying the rotation from az=el=0 local to fk5 boresight and 
    converts it into a rotation from az=el=0 local to galactic j2000 boresight coordinates
    '''

    if frame.type != core.G3FrameType.Scan:
        return
    gal_trans = G3VectorQuat()
    convert_ra_dec_trans_to_gal(frame[eq_trans_key], gal_trans)
    frame[out_key] = gal_trans


def AddAzElTrans(frame, az_key, el_key, out_key):
    '''
    Creates the transform for boresight pointing for az el based maps.
    
    Right now it's using the coordinate system where delta = -el because of implementation
    details.  All of the maps generated with this will be upside down in el.
    '''

    if frame.type != core.G3FrameType.Scan:
        return
    local_trans = G3VectorQuat()
    create_det_az_el_trans(frame[az_key], frame[el_key], local_trans)
    frame[out_key] = local_trans
