-----------
calibration
-----------

The calibration project contains data classes and analysis code for storing the following things:

*Physical properties of the bolometers*
	This includes the relative pointing offsets of the detectors on the focal plane, their polarization angles and efficiences, their bands, and the fabrication name (physical name) of the detector to which the channel corresponds. It does *not* include tuning-dependent parameters like time constants or responsivity. This information is stored in a ``BolometerPropertiesMap``, indexed by logical bolometer ID, in the Calibration frame. The standard name for this object is "BolometerProperties" and it is created from a composite of other calibration information.

