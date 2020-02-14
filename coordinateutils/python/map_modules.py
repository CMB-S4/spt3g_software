from spt3g import core
from spt3g.coordinateutils import G3SkyMapWeights, WeightType
import numpy as np

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

    frame['Q'] = np.zeros_like(frame['T'])
    frame['U'] = np.zeros_like(frame['T'])

    wgt_out = G3SkyMapWeights(frame['T'], weight_type=WeightType.Wpol)
    wgt_out.TT = wgt
    wgt_out.TQ = np.zeros_like(wgt)
    wgt_out.TU = np.zeros_like(wgt)
    wgt_out.QQ = np.ones_like(wgt)
    wgt_out.QU = np.zeros_like(wgt)
    wgt_out.UU = np.ones_like(wgt)

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

    wgt_out = G3SkyMapWeights(frame['T'], weight_type=WeightType.Wunpol)
    wgt_out.TT = wgt

    frame['Wunpol'] = wgt_out
