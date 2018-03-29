from spt3g import core
from spt3g.mapmaker import set_stokes_coupling, Mat1x3
import numpy as np

@core.cache_frame_data(type = core.G3FrameType.Map, bolo_props = 'BolometerProperties')
def ConvertTMapsToPolarized(frame, bolo_props = None):
    '''
    Converts individual weighted unpolarized maps to polarized versions of the same map.
    '''
    if frame.type != core.G3FrameType.Map:
        return
    assert(frame['Id'] in bolo_props)
    pol_angle = bolo_props[frame['Id']].pol_angle
    pol_efficiency = bolo_props[frame['Id']].pol_efficiency
    m = Mat1x3()
    set_stokes_coupling(pol_angle, pol_efficiency, 1, m)
    
    wgt = frame['Wunpol'].TT
    del frame['Wunpol']
    
    mp = frame['T'] * wgt
    mp.is_weighted = True
    del frame['T']
    
    frame['T'] = m[0] * mp
    frame['Q'] = m[1] * mp
    frame['U'] = m[2] * mp
    
    wgt_out = core.G3SkyMapWeights(mp, weight_type= core.WeightType.Wpol)
    wgt_out.TT = wgt * m[0] * m[0]
    wgt_out.TQ = wgt * m[0] * m[1]
    wgt_out.TU = wgt * m[0] * m[2]
    wgt_out.QQ = wgt * m[1] * m[1]
    wgt_out.QU = wgt * m[1] * m[2]
    wgt_out.UU = wgt * m[2] * m[2]

    frame['Wpol'] = wgt_out
    
@core.indexmod
def ConvertPolarizedMapsToT(frame):
    '''
    Converts individual weighted polarized maps to temperature-only versions of the same map.
    '''
    if frame.type != core.G3FrameType.Map:
        return

    wgt = frame['Wpol'].TT
    del frame['Wpol']
    del frame['Q']
    del frame['U']

    wgt_out = core.G3SkyMapWeights(frame['T'], weight_type= core.WeightType.Wunpol)
    wgt_out.TT = wgt

    frame['Wunpol'] = wgt_out
