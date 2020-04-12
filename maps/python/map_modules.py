from spt3g import core
from spt3g.maps import G3SkyMapWeights, get_mask_map
import numpy as np

@core.indexmod
def CompressMaps(frame, zero_nans=False):
    '''
    Compress all maps in a frame to their default sparse representation.
    Optionally remove NaN values as well.
    '''
    if frame.type != core.G3FrameType.Map:
        return

    for s in ['T', 'Q', 'U']:
        if s in frame:
            m = frame.pop(s)
            m.compress(zero_nans=zero_nans)
            frame[s] = m
    if 'Wunpol' in frame:
        m = frame.pop('Wunpol')
        m.TT.compress(zero_nans=zero_nans)
        frame['Wunpol'] = m
    if 'Wpol' in frame:
        m = frame.pop('Wpol')
        for s in ['TT', 'QQ', 'UU', 'TQ', 'TU', 'QU']:
            mm = getattr(m, s)
            mm.compress(zero_nans=zero_nans)
        frame['Wpol'] = m

@core.indexmod
def RemoveWeights(frame, zero_nans=False):
    '''
    Remove weights from input maps.  If zero_nans is `True`, empty pixles are
    skipped and pixels with zero weight are set to 0 instead of NaN.
    '''
    if frame.type != core.G3FrameType.Map:
        return

    if 'Wpol' not in frame and 'Wunpol' not in frame:
        return

    tmap = frame.pop('T')

    if 'Wpol' in frame:
        wmap = frame['Wpol']
        qmap = frame.pop('Q')
        umap = frame.pop('U')
    else:
        wmap = frame['Wunpol']
        qmap = None
        umap = None

    remove_weights(tmap, qmap, umap, wmap, zero_nans=zero_nans)

    frame['T'] = tmap
    if 'Wpol' in frame:
        frame['Q'] = qmap
        frame['U'] = umap

@core.indexmod
def ApplyWeights(frame):
    if frame.type != core.G3FrameType.Map:
        return

    if 'Wpol' not in frame and 'Wunpol' not in frame:
        return

    tmap = frame.pop('T')

    if 'Wpol' in frame:
        wmap = frame['Wpol']
        qmap = frame.pop('Q')
        umap = frame.pop('U')
    else:
        wmap = frame['Wunpol']
        qmap = None
        umap = None

    apply_weights(tmap, qmap, umap, wmap)

    frame['T'] = tmap
    if 'Wpol' in frame:
        frame['Q'] = qmap
        frame['U'] = umap

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
