from spt3g.calibration import BolometerProperties
from spt3g import core
import numpy as np
import re

__all__ = [
    "SplitByProperty",
    "SplitByBand",
    "SplitTimestreamsByBand",
    "SplitByWafer",
    "SplitByPixelType",
    "BandFormat",
    "get_band_units",
    "set_band_format",
    "band_to_string",
    "band_to_value",
    "extract_band_string",
    "extract_band_value",
]


# framework for handling frequency band formatting

class BandFormat:
    """
    Class for converting a frequency band between a quantity in G3Units and its
    string representation.

    Arguments
    ---------
    precision : int
        Float formatting precision for the band quantity.  If <=0, output will
        be an integer.  If >0, output will be a floating point number with this
        many decimal places of precision.
    units : str
        String name of the G3Units quantity in which to represent the frequency
        band, e.g. "GHz" or "MHz".  Must correspond to a valid attribute of the
        core.G3Units namespace.
    """

    def __init__(self, precision=0, units="GHz"):
        self.set_format(precision, units)

    def set_format(self, precision=0, units="GHz"):
        """
        Set the band format precision and units for converting between a
        quantity in G3Units and its string representation.

        Arguments
        ---------
        precision : int
            Float formatting precision for the band quantity.  If <=0, output
            will be an integer.  If >0, output will be a floating point number
            with this many decimal places of precision.
        units : str
            String name of the G3Units quantity in which to represent the
            frequency band, e.g. "GHz" or "MHz".  Must correspond to a valid
            attribute of the core.G3Units namespace.
        """
        self._precision = int(precision)
        assert hasattr(core.G3Units, units), "Invalid units {}".format(units)
        self._uvalue = getattr(core.G3Units, units)
        self._vformat = "%%.%df" % (precision if precision > 0 else 0)
        self._uformat = units
        self._format = "%s%s" % (self._vformat, units)
        prx = r"\.[0-9]{%d}" % precision if precision > 0 else ""
        self._pattern = "([0-9]+%s)%s" % (prx, units)
        self._wregex = re.compile("^" + self._pattern + "$")
        self._regex = re.compile(self._pattern)

    @property
    def units(self):
        """Units string used for formatting"""
        return self._uformat

    def to_string(self, value, include_units=True):
        """Convert a band value in G3Units to its string representation, using
        the appropriate precision and units name."""
        if not np.isfinite(value) or value < 0:
            return ""
        value = np.round(value / self._uvalue, self._precision)
        if not include_units:
            return self._vformat % value
        return self._format % value

    def to_value(self, string):
        """Convert a band string to a value in G3Units, or raise a ValueError
        if no match is found."""
        m = self._wregex.match(string)
        if not m:
            raise ValueError("Invalid band {}".format(string))
        return float(m.group(1)) * self._uvalue

    def extract_string(self, string, include_units=True):
        """Return the band substring from the input string, or None if not
        found."""
        m = self._regex.match(string)
        if not m:
            return None
        if not include_units:
            return m.group(1)
        return m.group(0)

    def extract_value(self, value):
        """Return the band in G3Units extracted from the input string, or None
        if not found."""
        s = extract_string(value, include_units=False)
        if not s:
            return None
        return float(v) * self._uvalue


# global instance used by functions below
_band_format = BandFormat()


@core.usefulfunc
def get_band_units():
    """Return the units string used for formatting frequency band values."""
    return _band_format.units


@core.usefulfunc
def set_band_format(precision, units):
    _band_format.set_format(precision, units)
set_band_format.__doc__ = BandFormat.set_format.__doc__


@core.usefulfunc
def band_to_string(value, include_units=True):
    return _band_format.to_string(value, include_units=include_units)
band_to_string.__doc__ = BandFormat.to_string.__doc__


@core.usefulfunc
def band_to_value(string):
    return _band_format.to_value(string)
band_to_value.__doc__ = BandFormat.to_value.__doc__


@core.usefulfunc
def extract_band_string(string, include_units=True):
    return _band_format.extract_string(string, include_units=include_units)
extract_band_string.__doc__ = BandFormat.extract_string.__doc__


@core.usefulfunc
def extract_band_value(string):
    return _band_format.extract_value(string)
extract_band_value.__doc__ = BandFormat.extract_value.__doc__


# monkeypatch useful property attributes
def band_string(self):
    """String representation of frequency band center"""
    return _band_format.to_string(self.band)
BolometerProperties.band_string = property(band_string)

def band_vstring(self):
    """String representation of frequency band center, without units name"""
    return _band_format.to_string(self.band, include_units=False)
BolometerProperties.band_vstring = property(band_vstring)


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
