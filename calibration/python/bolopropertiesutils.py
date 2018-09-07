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
class SplitByProperty(object):
    '''
    Take an input G3FrameObject-derivative Map keyed by bolometer name and
    split it into several based on the property of the detectors as given by
    the BolometerProperties key.
    Return the same type of maps as the one it was handed, e.g.
    G3TimestreamMap, G3MapInt, etc.
    '''
    def __init__(self, input='CalTimestreams', property=None,
                 converter=lambda x: str(x), property_list=None,
                 output_root=None, bpm='BolometerProperties'):
        '''
        Split the input map given by input into several output
        maps named output_root + key (e.g. CalTimestreams + converter(property)) with
        the default options).

        Arguments
        ---------
        input : str
            Key name of the input map to split.
        property : str
            Attribute name to extract from the BolometerProperties object.
            Required.
        converter : callable
            Callable (function) for converting the property to its corresponding
            string name.  Returns a string representation of the input argument,
            or None if the argument is invalid.
        property_list : list of properties
            Properties to include in the output keys.  Entries that are not strings
            will be converted to strings using `converter`.
            If property_list is not None, use only the names in the
            list (possibly writing empty timestream maps to the frame). Otherwise,
            creates maps for every that exists in the input.
        output_root : str
            Prefix for the output keys.
            If None (default), use `input` as the output root.
        bpm : str
            The key name of the BolometerPropertiesMap from which to extract
            the requested `property` for splitting the input map.
        '''
        if property is None:
            core.log_fatal("Property is a required argument")
        self.bpmattr = property

        if converter is None:
            core.log_fatal("Converter is a required argument")
        self.converter = converter

        self.input = input
        self.output_root = output_root if output_root is not None else input
        if property_list is not None:
            self.props = [converter(x) if not isinstance(x, str) else x
                          for x in property_list]
        else:
            self.props = None
        self.bpmkey = bpm
        self.bpm = None

    def __call__(self, frame):
        if self.bpmkey in frame:
            self.bpm = frame[self.bpmkey]

        if self.input not in frame:
            return

        inmap = frame[self.input]
        out = {}
        if self.props is not None:
            for prop in self.props:
                out[prop] = type(inmap)()

        for b in inmap.keys():
            try:
                prop = self.converter(getattr(self.bpm[b], self.bpmattr))
            except KeyError:
                continue
            if prop not in out:
                if self.props is None and prop is not None:
                    out[prop] = type(inmap)()
                else:
                    continue
            out[prop][b] = inmap[b]

        for prop in out.keys():
            frame['%s%s' % (self.output_root, prop)] = out[prop]


@core.indexmod
class SplitByBand(SplitByProperty):
    '''
    Take an input G3FrameObject-derivative Map keyed by bolometer name and
    split it into several based on the bands of the detectors as given by
    the BolometerProperties key.
    Return the same type of maps as the one it was handed, e.g.
    G3TimestreamMap, G3MapInt, etc.
    '''
    def __init__(self, input='CalTimestreams', output_root=None,
                 bands=None, bpm='BolometerProperties'):
        '''
        Split the input map given by input into several output
        maps named output_root + band + GHz (e.g. CalTimestreams150GHz with
        the default options). If bands is not None, use only the bands in the
        list (possibly writing empty timestream maps to the frame). Otherwise,
        creates maps for every band that exists in the input. Setting bpm
        to a non-default value causes this to get its band mapping from an
        alternative data source.
        '''
        def converter(band):
            if isinstance(band, str):
                return band
            if math.isnan(band) or math.isinf(band):
                return None
            if band < 0:
                return None
            return '%dGHz' % int(band/core.G3Units.GHz)

        super(SplitByBand, self).__init__(
            input=input, output_root=output_root, property_list=bands,
            bpm=bpm, property='band', converter=converter)


@core.indexmod
class SplitByWafer(SplitByProperty):
    '''
    Take an input G3FrameObject-derivative Map keyed by bolometer name and
    split it into several based on the wafers of the detectors as given by
    the BolometerProperties key.
    Return the same type of maps as the one it was handed, e.g.
    G3TimestreamMap, G3MapInt, etc.
    '''
    def __init__(self, input='CalTimestreams', output_root=None,
                 wafers=None, bpm='BolometerProperties'):
        '''
        Split the input map given by input into several output
        maps named output_root + wafer (e.g. CalTimestreamsW172 with
        the default options). If wafers is not None, use only the wafers in the
        list (possibly writing empty timestream maps to the frame). Otherwise,
        creates maps for every wafer that exists in the input. Setting bpm
        to a non-default value causes this to get its wafer mapping from an
        alternative data source.
        '''
        def converter(wafer):
            if wafer is None:
                return None
            return str(wafer).capitalize()

        super(SplitByWafer, self).__init__(
            input=input, output_root=output_root, property_list=wafers,
            bpm=bpm, property='wafer_id', converter=converter)
