from spt3g.calibration import BolometerProperties
from spt3g import core
import math

__all__ = ['SplitByProperty', 'SplitByBand', 'SplitTimestreamsByBand',
           'SplitByWafer', 'SplitByPixelType']

@core.indexmod
class SplitByProperty(object):
    '''
    Take an input G3FrameObject-derivative Map keyed by bolometer name and
    split it into several based on the property of the detectors as given by
    the BolometerProperties key.
    Return the same type of maps as the one it was handed, e.g.
    G3TimestreamMap, G3MapInt, etc.
    '''
    def __init__(self, input='CalTimestreams', property=None, property_list=None,
                 output_root=None, bpm='BolometerProperties'):
        '''
        Split the input map given by input into several output
        maps named output_root + key (e.g. CalTimestreams + str(property)) with
        the default options).

        Arguments
        ---------
        input : str
            Key name of the input map to split.
        property : str
            Attribute name to extract from the BolometerProperties object.
            Required.
        property_list : list of properties
            Properties to include in the output keys.  Entries that are not strings
            will be converted to strings using the `SplitByProperty.converter` method.
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

        self.input = input
        self.output_root = output_root if output_root is not None else input
        if property_list is not None:
            self.props = [self.converter(x) if not isinstance(x, str) else x
                          for x in property_list]
        else:
            self.props = None
        self.bpmkey = bpm
        self.bpm = None

    @staticmethod
    def converter(prop):
        """
        Function for converting the property to its corresponding string name.
        Returns a string representation of the input argument, or None if the
        argument is invalid.

        Overload this function in subclasses of SplitByProperty to change
        how attributes are parsed into their string representations.
        """
        if prop is None:
            return None
        return str(prop)

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
        super(SplitByBand, self).__init__(
            input=input, output_root=output_root, property_list=bands,
            bpm=bpm, property='band')

    @staticmethod
    def converter(band):
        if isinstance(band, str):
            return band
        if math.isnan(band) or math.isinf(band):
            return None
        if band < 0:
            return None
        return '%dGHz' % int(band/core.G3Units.GHz)


@core.indexmod
class SplitTimestreamsByBand(SplitByBand):

    def __init__(self, input='CalTimestreams', output_root=None,
                 bands=None, bpm='BolometerProperties'):
        core.log_warn("SplitTimestreamsByBand is deprecated, use SplitByBand instead")
        super(SplitTimestreamsByBand, self).__init__(
            input=input, output_root=output_root, bands=bands, bpm=bpm)


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
        super(SplitByWafer, self).__init__(
            input=input, output_root=output_root, property_list=wafers,
            bpm=bpm, property='wafer_id')

    @staticmethod
    def converter(wafer):
        if wafer is None:
            return None
        return str(wafer).capitalize()


@core.indexmod
class SplitByPixelType(SplitByProperty):
    '''
    Take an input G3FrameObject-derivative Map keyed by bolometer name and
    split it into several based on the pixel types of the detectors as given by
    the BolometerProperties key.
    Return the same type of maps as the one it was handed, e.g.
    G3TimestreamMap, G3MapInt, etc.
    '''
    def __init__(self, input='CalTimestreams', output_root=None,
                 types=None, bpm='BolometerProperties'):
        '''
        Split the input map given by input into several output
        maps named output_root + wafer (e.g. CalTimestreamsW172 with
        the default options). If wafers is not None, use only the wafers in the
        list (possibly writing empty timestream maps to the frame). Otherwise,
        creates maps for every wafer that exists in the input. Setting bpm
        to a non-default value causes this to get its wafer mapping from an
        alternative data source.
        '''
        super(SplitByPixelType, self).__init__(
            input=input, output_root=output_root, property_list=types,
            bpm=bpm, property='pixel_type')

    @staticmethod
    def converter(pixel_type):
        if pixel_type is None:
            return None
        if not pixel_type:
            return None
        pixel_type = str(pixel_type)
        if pixel_type.lower() == 'n/a':
            return None
        if pixel_type.islower():
            return pixel_type.capitalize()
        return pixel_type
