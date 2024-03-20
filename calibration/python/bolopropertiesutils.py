from spt3g.calibration import BolometerProperties
from spt3g import core
import numpy as np
import re

__all__ = ['SplitByProperty', 'SplitByBand', 'SplitTimestreamsByBand',
           'SplitByWafer', 'SplitByPixelType']


# framework for handling frequency band formatting

BolometerProperties._band_precision = 0
BolometerProperties._band_units = "GHz"

def get_band_precision(cls):
    """Bolometer frequency band precision for string formatting"""
    return cls._band_precision
def set_band_precision(cls, precision):
    cls._band_precision = int(precision)
BolometerProperties.band_precision = property(get_band_precision, set_band_precision)

def get_band_units(cls):
    """Bolometer units name for string formatting"""
    return cls._band_units
def set_band_units(cls, units):
    assert hasattr(core.G3Units, units), "Invalid units {}".format(units)
    cls._band_units = units
BolometerProperties.band_units = property(get_band_units, set_band_units)

def band_to_string(cls, band):
    """Convert a band number in G3Units to its string representation, using the
    appropriate precision and units name."""
    precision = cls.band_precision
    units = cls.band_units
    band = np.round(band / getattr(core.G3Units, units), precision)
    return "%%.%df%s" % (precision, units) % band
BolometerProperties.band_to_string = classmethod(band_to_string)

def band_string(self):
    """Band string representation"""
    if not np.isfinite(self.band) or self.band < 0:
        return ""
    return self.band_to_string(self.band)
BolometerProperties.band_string = property(band_string)

def band_to_num(cls, band):
    """Convert a band string to a number in G3Units, or raise a ValueError if no
    match is found."""
    m = cls.band_regex.match(band)
    if not m:
        raise ValueError("Invalid band {}".format(band))
    return float(m.group(1)) * getattr(core.G3Units, units)
BolometerProperties.band_to_num = classmethod(band_to_num)

def band_regex(cls):
    """Compiled regular expression for identifying a band string, whose first
    match group returns the numerical portion of the string."""
    precision = cls.band_precision
    units = cls.band_units
    prx = "\\.[0-9]{{}}".format(precision) if precision else ""
    return re.compile("([0-9]+{}){}".format(prx, units))
BolometerProperties.band_regex = property(band_regex)

def band_string_from_string(cls, s):
    """Return the band substring from the input string, or None if not found."""
    m = cls.band_regex.search(s)
    if not m:
        return None
    return m.group(0)
BolometerProperties.band_string_from_string = classmethod(band_string_from_string)

def band_from_string(cls, s):
    """Return the band in G3Units extracted from the input string, or None if
    not found."""
    m = cls.band_regex.search(s)
    if not m:
        return None
    return float(m.group(1)) * getattr(core.G3Units, cls.band_units)
BolometerProperties.band_from_string = classmethod(band_from_string)


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
                 output_root=None, bpm='BolometerProperties', drop_empty=False):
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
        drop_empty : bool
            If True, drop output maps that don't contain any bolometers.
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
        self.drop_empty = drop_empty

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

        for prop, outmap in out.items():
            if not len(outmap.keys()) and self.drop_empty:
                continue
            frame['%s%s' % (self.output_root, prop)] = outmap


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
                 bands=None, bpm='BolometerProperties', drop_empty=False):
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
            bpm=bpm, property='band_string', drop_empty=drop_empty)

    @staticmethod
    def converter(band_string):
        if band_string is None:
            return None
        if not band_string:
            return None
        return str(band_string)


@core.indexmod
class SplitTimestreamsByBand(SplitByBand):

    def __init__(self, input='CalTimestreams', output_root=None,
                 bands=None, bpm='BolometerProperties', drop_empty=False):
        core.log_warn("SplitTimestreamsByBand is deprecated, use SplitByBand instead")
        super(SplitTimestreamsByBand, self).__init__(
            input=input, output_root=output_root, bands=bands, bpm=bpm, drop_empty=drop_empty)


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
                 wafers=None, bpm='BolometerProperties', drop_empty=False):
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
            bpm=bpm, property='wafer_id', drop_empty=drop_empty)

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
                 types=None, bpm='BolometerProperties', drop_empty=False):
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
            bpm=bpm, property='pixel_type', drop_empty=drop_empty)

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
