-----------
calibration
-----------

The calibration project contains data classes and analysis code for storing the following things:

*Physical properties of the bolometers*
	This includes the relative pointing offsets of the detectors on the focal plane, their polarization angles and efficiences, their bands, and the fabrication name (physical name) of the detector to which the channel corresponds. It does *not* include tuning-dependent parameters like time constants or responsivity. This information is stored in a ``BolometerPropertiesMap``, indexed by logical bolometer ID, in the Calibration frame. The standard name for this object is "BolometerProperties" and it is created from a composite of other calibration information.

	The frequency band for a bolometer is often represented by a string name with units.  The ``BolometerProperties`` class provides some functions for converting numerical bands (in G3Units) to their string representation.  By default, the band string is formatted with precision 0 (i.e. an integer) in units of GHz; these can be modified by using the ``BolometerProperties.set_band_format()`` class method.  The ``.band_string`` property of a ``BolometerProperties`` instance returns the string representation of the given bolometer's frequency band.  The ``BolometerProperties.band_to_string`` and ``BolometerProperties.band_to_num`` class methods can be used to convert numerical bands to their string representation and vice versa.  Library code should take care to use these methods wherever possible.  The ``BolometerProperties.band_string_from_string`` and ``BolometerProperties.band_from_string`` class methods can be used to extract band information in string or numerical form from other strings.
