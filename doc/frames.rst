------
Frames
------

The basic unit of data processing is the frame. Each frame is a free-form mapping from strings to data of a type derived from G3FrameObject. In general, they behave as Python dictionaries. Wrappers are provided for Python plain-old-data types (numbers, strings, etc.) so that they can be used directly.

Frames have fast serialization to and from disk (see the G3Reader and G3Writer modules in :doc:`fileio`) with built-in data integrity protection. This serialization is accessible programmatically in Python through the Python pickle libraries, which are overridden to use the internal serialization.

Each frame has a type defined under the ``core.G3FrameType`` namespace. These are meant to indicate different types of data (calibration vs. scan data, for example) and, in general, many types of frames will be interleaved in the same data stream. A good general rule for whether data should be in a different frame type is to consider the rates at which the data change: data that change at the same speed (e.g. bolometer data and pointing) should share a frame, while data that change at different speeds (e.g. calibration constants and bolometer data) should be separated.

A brief description of the intention of each frame type follows along with a table containing a representative minimal set of data likely to be contained by a frame of a given type. The list is neither exhaustive nor truly minimal: other data can and will be present and some of the data listed here may have been removed or renamed. Neither is this an exhaustive list of frame types: any single character code can be [ab]used as a frame type for special purpose tools.

.. contents:: Contents

Scan
====

Scan frames are the most common unit of data for analysis. They contain the data from a single left-to-right (or right-to-left) scan of the telescope. Where quantities are bolometer-indexed, the names match those in the Wiring map.

===================	====================	===========
Key			Type			Description
===================	====================	===========
ACUStatus		ACUStatusVector		Sequences of gcp.ACUStatus objects containing the ACU state readout during the scan
SourceName		G3String		Name of the source being observed
BoresightAz		G3Timestream		Telescope azimuth as a function of time
BoresightEl		G3Timestream		Telescope elevation as a function of time
RawTimestreams_I	G3TimestreamMap		Raw (counts) timestreams from detectors, indexed by bolometer, I modulation
RawTimestreams_Q	G3TimestreamMap		Raw (counts) timestreams from detectors, indexed by bolometer, Q modulation
CalTimestreams		G3TimestreamMap		Bolometer timestreams with some useful calibration applied (not present in raw data)
TimestreamWeights	G3MapDouble		Floating point weights for each timestream, indexed by bolometer (not present in raw data)
ScanNumber		G3Int			Sequence number of the current scan within an observation
Flags			G3MapVectorString	List of detectors with exceptional conditions. Indexed by bolometer ID. Each element is a list of strings that provides a description of the condition[s] that caused the detector to be flagged (not present in raw data).
Turnaround		G3Bool			Present and set to True if scan is a turnaround or otherwise is not a constant velocity scan. Absent otherwise.
TrakerStatus		TrackerStatus		Composite tracker status information at all points during the scan.
DfMuxHousekeeping	DfMuxHousekeepingMap		Tree containing housekeeping data and configuration from the boards connected to the readout system, indexed by board IP address
CalibratorOn		G3Timestream		Timestream of the calibrator sync signal. Set to 1 if the sync line is high and 0 otherwise. NaN if no readout.
GCPFeatureBits		G3VectorString		Strings describing the GCP flags field.
===================	====================	===========

PipelineInfo
============

PipelineInfo frames contain information on the processing applied to make the remainder of the data. They typically occur at the beginning of the data stream and are automatically inserted by G3Pipeline, if not previously present. The frame will contain at least one element of type G3PipelineInfo that gives information on the version of the software that made the data and the configuration of all modules and/or segments added to the pipeline. Other data (G3PipelineInfo from previous scripts, configuration information added by other modules) may also be present.

==============================		===============	===========
Key					Type		Description
==============================		===============	===========
22-May-2019:18:42:15.335969000		G3PipelineInfo	Configuration information of the pipeline that produced an object. Added automatically. The key name is a timestamp indicating the time at which the software began running.
==============================		===============	===========

Timepoint
=========

Timepoint frames include sample-by-sample data from the mux system. They are emitted directly by the DfMuxBuilder as part of primary data acquisition and are generally not seen outside the South Pole or the lab. These are transformed into Scan_ frames by DfMuxCollator.

===================	===============	===========
Key			Type		Description
===================	===============	===========
EventHeader		G3Time		Time of the sample
DfMux			DfMuxMetaSample	Tree containing data from the boards connected to the readout system, indexed by board IP address
===================	===============	===========

Housekeeping
============

Contains housekeeping data. Issued periodically when housekeeping data is taken. Like Timepoint_ frames, these are rolled into Scan_ frames during processing and do not appear in general in stored data.

===================	====================		===========
Key			Type				Description
===================	====================		===========
DfMuxHousekeeping	DfMuxHousekeepingMap		Tree containing housekeeping data and configuration from the boards connected to the readout system, indexed by board IP address
===================	====================		===========

Map
===

Contains either the result of the map maker or the input to simulation.

==========================	======================	===========
Key				Type			Description
==========================	======================	===========
Id				G3String		A string identifying the map for the various processing steps
T				G3SkyMap		A map storing the intensity information (could be sky intensity or sky intensity x weight).  In the case of maps that store abstract information like apodization masks or point source masks, the data will also be stored under the T key.  The motivation being it makes it easy to have G3Modules operating on maps also work on these.
Q				G3SkyMap		A map storing the Q polarization information (could be sky Q or sky Q x weight)
U				G3SkyMap		A map storing the u polarization information (could be sky U or sky U x weight)
Wpol				G3SkyMapWeights		If the frame contains polarized information, this stores the t/q/u covariances scaled by the individual detector weights
Wunpol				G3SkyMapWeights		This stores the unpolarized weight information
==========================	======================	===========



Calibration
===========

This frame contains all measured calibration information (pointing, response, etc.) that may change when remeasured. It does *not* include static information describing how the instrument is set up (see Wiring_ below).

==========================	======================	===========
Key				Type			Description
==========================	======================	===========
BolometerProperties		BolometerPropertiesMap	Measured non-configuration-dependent calibration properties of the instrument (pointing, pol efficiency, etc.), indexed by bolometer
NominalBolometerProperties	BolometerPropertiesMap	As above, but what those properties were meant to be.
TimeConst			G3MapDouble		Time constants of the detectors. These can change with the bias point. Should perhaps be moved to the InstrumentStatus frame.
RCW38FluxCalibration		G3MapDouble		Observed flux of RCW38 per detector as a fraction of the calibrator response.
CalibratorResponse		G3MapDouble		Observed response to the most recent calibrator observation for each detector in Watts.
CalibratorResponseSN		G3MapDouble		Signal to noise (in sigma) of the most recent calibrator observation in sigma.
==========================	======================	===========

Observation
===========

Indicates global observation parameters. Changes at the beginning of a new observation, though, as with all metadata, repeat observation frames may appear during processing.

=========================	======================	===========
Key				Type			Description
=========================	======================	===========
SourceName			G3String		Name of the source being observed
ObservationNumber		G3Int			Sequence number of the current observation since we started recording such things
=========================	======================	===========

Wiring
======

Gives the description of how the system is wired: notably, the connection between board serial number, module, channel and a bolometer ID.

=============	==============	======================================
Key		Type		Description
=============	==============	======================================
WiringMap	DfMuxWiringMap	Digest of the pydfmux channel mappings
ReadoutSystem	G3String	Description of the type of readout system employed. Set to "DfMux" for SPTpol-style readout and "ICE" for 3G-style readout.
=============	==============	======================================

GcpSlow
=======

Holds all the GCP data sampled once per second. Like Timepoint_ and Housekeeping_ frames, these data are consolidated in the Scan frames and do not appear in final data products. The data stored here begins as a strict copy of the GCP register map (see the GCP documentation for details, a few notable entries are summarized below). A few other keys are added transiently in the course of generating Scan_ frames.

=========	================	======================================
Key		Type			Description
=========	================	======================================
array		G3MapFrameObject	Most of the data stored by GCP
antenna0	G3MapFrameObject	Telescope pointing information
=========	================	======================================



EndProcessing
=============

EndProcessing is a special-purpose frame emitted implicitly by G3Pipeline at the end of processing. No further frames will occur after this and reception of an EndProcessing frame is intended as a signal to modules that they should clear any internally buffered data and clean up.

EndProcessing frames should, in general, contain no data.

