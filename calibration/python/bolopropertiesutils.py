from spt3g.calibration import BolometerProperties
from spt3g import core
import math

def InvertBoloProperties(frame, 
                         bolo_props_key = 'BolometerProperties',
                         wafer_map_key = 'WaferMap',
                         squid_map_key = 'SquidMap',
                         band_map_key = 'BandMap'):
    if frame.type != core.G3FrameType.Calibration:
        return

    def _add_myself(k, val, mp):
        if not val in mp:
            mp[val] = []
        mp[val] = k

    bolo_props = frame[bolo_props_key]
    
    wafer_map = core.G3MapVectorString()
    squid_map = core.G3MapVectorString()
    band_map = core.G3MapVectorString()

    for k in bolo_props.keys():
        _add_myself(k, bolo_props[k].wafer_id, wafer_map)
        _add_myself(k, bolo_props[k].squid_id, squid_map)
        _add_myself(k, bolo_props[k].band, band_map)

    frame[wafer_map_key] = wafer_map
    frame[squid_map_key] = squid_map
    frame[band_map_key] = band_map
    
@core.indexmod
class SplitTimestreamsByBand(object):
    '''
    Take an input G3TimestreamMap and split it into several based on
    the bands of the detectors as given by the BolometerProperties
    key.
    '''
    def __init__(self, input='CalTimestreams', output_root='CalTimestreams',
      bands=None, bpm='BolometerProperties'):
        '''
        Split the input timestream map given by input into several output
        maps named output_root + band + GHz (e.g. CalTimestreams150GHz with
        the default options). If bands is not None, use only the bands in the
        list (possibly writing empty timestream maps to the frame). Otherwise,
        creates maps for every band that exists in the input. Setting bpm
        to a non-default value causes this to get its band mapping from an
        alternative data source.
        '''
        self.input = input
        self.output_root = output_root
        self.bands = bands
        self.bpmkey = bpm
        self.bpm = None

    def __call__(self, frame):
        if self.bpmkey in frame:
            self.bpm = frame[self.bpmkey]

        if self.input in frame:
            out = {}
            if self.bands is not None:
                for band in self.bands:
                    out[band] = core.G3TimestreamMap()
            
            for b, ts in frame[self.input].iteritems():
                try:
                    band = self.bpm[b].band
                except KeyError:
                    continue
                if band not in out:
                    if self.bands is None and \
                      not (math.isnan(band) or math.isinf(band)):
                        out[band] = core.G3TimestreamMap()
                    else:
                        continue
                out[band][b] = ts

            for band in out.keys():
                frame['%s%dGHz' % (self.output_root,
                  int(band/core.G3Units.GHz))] = out[band]
    
