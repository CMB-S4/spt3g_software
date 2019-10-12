from spt3g import core
from spt3g.maps import G3SkyMapWeights
import numpy as np

@core.indexmod
def ConvertTMapsToPolarized(frame):
    '''
    Converts individual unpolarized maps to polarized versions of the same map.

    This module is only a shim that creates null Q and U maps and populates
    a properly invertible Wpol array from the TT Wunpol weights.
    '''
    if frame.type != core.G3FrameType.Map:
        return

    wgt = frame['Wunpol'].TT
    del frame['Wunpol']

    frame['Q'] = 0 * frame['T']
    frame['U'] = 0 * frame['T']

    wgt_out = G3SkyMapWeights(frame['T'], ispolarized=True)
    wgt_out.TT = wgt
    wgt_out.TQ = 0 * wgt
    wgt_out.TU = 0 * wgt
    wgt_out.QQ = 0 * wgt + 1.0
    wgt_out.QU = 0 * wgt
    wgt_out.UU = 0 * wgt + 1.0

    frame['Wpol'] = wgt_out

@core.indexmod
def ConvertPolarizedMapsToT(frame):
    '''
    Converts individual polarized maps to temperature-only versions of the same map.
    '''
    if frame.type != core.G3FrameType.Map:
        return

    wgt = frame['Wpol'].TT
    del frame['Wpol']
    del frame['Q']
    del frame['U']

    wgt_out = G3SkyMapWeights(frame['T'], ispolarized=False)
    wgt_out.TT = wgt

    frame['Wunpol'] = wgt_out
