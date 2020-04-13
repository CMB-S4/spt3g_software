from spt3g import core, maps
import numpy as np

__all__ = [
    'CompressMaps',
    'RemoveWeights',
    'ApplyWeights',
    'FlattenPol',
    'ConvertTMapsToPolarized',
    'ConvertPolarizedMapsToT',
]

@core.indexmod
def CompressMaps(frame, zero_nans=False):
    '''
    Compress all maps in a frame to their default sparse representation.
    Optionally remove NaN values as well.
    '''
    for s in ['T', 'Q', 'U', 'Wunpol', 'Wpol']:
        if s in frame:
            m = frame.pop(s)
            m.compress(zero_nans=zero_nans)
            frame[s] = m

@core.indexmod
def RemoveWeights(frame, zero_nans=False):
    '''
    Remove weights from input maps.  If zero_nans is `True`, empty pixles are
    skipped and pixels with zero weight are set to 0 instead of NaN.  Operation
    is performed in place to minimuze memory use.
    '''
    if 'Wpol' not in frame and 'Wunpol' not in frame:
        return

    tmap = frame.pop('T')

    if 'Wpol' in frame:
        wmap = frame['Wpol']
        qmap = frame.pop('Q')
        umap = frame.pop('U')
        maps.remove_weights(tmap, qmap, umap, wmap, zero_nans=zero_nans)
    else:
        wmap = frame['Wunpol']
        maps.remove_weights_t(tmap, wmap, zero_nans=zero_nans)

    frame['T'] = tmap
    if 'Wpol' in frame:
        frame['Q'] = qmap
        frame['U'] = umap

@core.indexmod
def ApplyWeights(frame):
    '''
    Apply weights to the input maps.  The operation is performed in place to
    minimize memory use.
    '''
    if 'Wpol' not in frame and 'Wunpol' not in frame:
        return

    tmap = frame.pop('T')

    if 'Wpol' in frame:
        wmap = frame['Wpol']
        qmap = frame.pop('Q')
        umap = frame.pop('U')
        maps.apply_weights(tmap, qmap, umap, wmap)
    else:
        wmap = frame['Wunpol']
        maps.apply_weights_t(tmap, wmap)

    frame['T'] = tmap
    if 'Wpol' in frame:
        frame['Q'] = qmap
        frame['U'] = umap

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

    if 'Q' not in frame or 'U' not in frame:
        return

    qmap, umap = frame['Q'], frame['U']
    if not isinstance(qmap, maps.FlatSkyMap) or not isinstance(umap, maps.FlatSkyMap):
        return

    maps.flatten_pol(qmap, umap, invert=invert)
    return frame

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
    mask = maps.get_mask_map(wgt)

    wgt_out = maps.G3SkyMapWeights(frame['T'], polarized=True)
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

    wgt_out = maps.G3SkyMapWeights(frame['T'], polarized=False)
    wgt_out.TT = wgt

    frame['Wunpol'] = wgt_out
