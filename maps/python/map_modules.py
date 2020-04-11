from spt3g import core
from spt3g.maps import G3SkyMapWeights, get_mask_map, zero_map_nans
import numpy as np

@core.indexmod
def MakeMapsSparse(frame):
    '''
    Makes all maps in a frame sparse.
    '''
    if frame.type != core.G3FrameType.Map:
        return

    for s in ['T', 'Q', 'U']:
        if s in frame:
            m = frame.pop(s)
            m.sparse = True
            frame[s] = m
    if 'Wunpol' in frame:
        m = frame.pop('Wunpol')
        m.TT.sparse = True
        frame['Wunpol'] = m
    if 'Wpol' in frame:
        m = frame.pop('Wpol')
        for s in ['TT', 'QQ', 'UU', 'TQ', 'TU', 'QU']:
            getattr(m, s).sparse = True
        frame['Wpol'] = m

@core.indexmod
def ZeroMapNans(frame, check_weights=False):
    '''
    Convert NaN values in the input maps and (optionally) weights to zeros.
    '''
    if frame.type != core.G3FrameType.Map:
        return

    for s in ['T', 'Q', 'U']:
        if s in frame:
            m = frame.pop(s)
            frame[s] = zero_map_nans(m)

    if not check_weights:
        return

    if 'Wunpol' in frame:
        m = frame.pop('Wunpol')
        m.TT = zero_map_nans(m.TT)
        frame['Wunpol'] = m

    if 'Wpol' in frame:
        m = frame.pop('Wpol')
        for s in ['TT', 'QQ', 'UU', 'TQ', 'TU', 'QU']:
            setattr(m, s, zero_map_nans(getattr(m, s)))
        frame['Wpol'] = m

@core.indexmod
def ConvertTMapsToPolarized(frame):
    '''
    Converts individual unpolarized maps to polarized versions of the same map.

    This module is only a shim that creates null Q and U maps and populates
    a properly invertible Wpol array from the TT Wunpol weights.
    '''
    if frame.type != core.G3FrameType.Map or 'Wunpol' not in frame:
        return

    wgt = frame['Wunpol'].TT
    del frame['Wunpol']

    frame['Q'] = frame['T'].Clone(False)
    frame['U'] = frame['T'].Clone(False)
    mask = get_mask_map(wgt)

    wgt_out = G3SkyMapWeights(frame['T'], polarized=True)
    wgt_out.TT = wgt
    wgt_out.TQ = wgt.Clone(False)
    wgt_out.TU = wgt.Clone(False)
    wgt_out.QQ = mask
    wgt_out.QU = wgt.Clone(False)
    wgt_out.UU = mask.Clone(True)

    frame['Wpol'] = wgt_out

@core.indexmod
def ConvertPolarizedMapsToT(frame):
    '''
    Converts individual polarized maps to temperature-only versions of the same map.
    '''
    if frame.type != core.G3FrameType.Map or 'Wpol' not in frame:
        return

    wgt = frame['Wpol'].TT
    del frame['Wpol']
    del frame['Q']
    del frame['U']

    wgt_out = G3SkyMapWeights(frame['T'], polarized=False)
    wgt_out.TT = wgt

    frame['Wunpol'] = wgt_out
