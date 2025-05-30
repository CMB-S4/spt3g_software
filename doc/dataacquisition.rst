----------------
Data Acquisition
----------------

This software provides, in addition to data processing, facilities for data acquisition. The goal is to have a unified framework for simulations, data acquisition, and analysis.

.. contents:: Contents
   :local:

DfMux
=====

Bolometer Data
~~~~~~~~~~~~~~

The :doc:`dfmux <moddoc_dfmux>` project provides two classes that can be used for event-driven acquisition of bolometer timestreams from IceBoards. The core class, :py:class:`~spt3g.dfmux.DfMuxCollector`, listens for data from one or more IceBoards, passing the resulting data to a :py:class:`~spt3g.dfmux.DfMuxBuilder` object to which the collected data is reported. The :py:class:`~spt3g.dfmux.DfMuxBuilder` is a subclass of :py:class:`~spt3g.core.G3EventBuilder` (see G3EventBuilder_ below) that assembles mux packets from one or more :py:class:`~spt3g.dfmux.DfMuxCollector` objects into frames, with one sample from each readout channel per frame. To do this, it is passed either the number of boards to expect for a complete sample or a list of the serial numbers of boards to use.

.. code-block:: python

	pipe = core.G3Pipeline()
	builder = dfmux.DfMuxBuilder([136])
	collector = dfmux.DfMuxCollector(builder, ["iceboard0136.local"])
	collector.Start()
	pipe.Add(builder)

The above example listens for data from the IceBoard with serial number 136 using SCTP (see below), passing the resulting data to the specified DfMuxBuilder. A list of serial numbers in a given pydfmux hardware map can be obtained using a command like this:

.. code-block:: python

        serials = [int(board.serial) for board in hwm.query(pydfmux.IceBoard)]


This library supports data acquisition using either of two data transport mechanisms: SCTP (newer firmware than 11.3 only) and multicast UDP.

SCTP
____

Extremely new (as-yet unreleased) IceBoard firmwares support using SCTP as a data transport. In this mode, a :py:class:`~spt3g.dfmux.DfMuxCollector` object connects to the IceBoard, opening a data connection over which streaming samples are transported to the DAQ computer. This connection is error-tolerant and point-to-point, so the collector must be passed a list of boards to listen to. SCTP is used when the collector is created with a list of board hostnames and the builder argument *first*:

.. code-block:: python

	collector = dfmux.DfMuxCollector(builder, ["iceboard0136.local"])

Note that, when using SCTP for data transport on Linux, you may need to load the sctp kernel module by running ``modprobe sctp``. On versions of the Linux kernel earlier than 4.16, you may also need to instantiate one DfMuxCollector per board. (This is the default behavior of ``record_bolodata.py``).


UDP
___

The default mode of data acquisition is multicast UDP, which makes the data acquisition system passive. When using multicast UDP (mandatory on firmwares 11.3 and earlier, otherwise optional), the :py:class:`~spt3g.dfmux.DfMuxCollector` must be passed the IP address of an interface *on the DAQ computer* on which to listen for detector data and, optionally, a list of board serial numbers. (The list of board serial numbers is optional when using only one Ethernet interface, but *must* be passed if using Linux on a system with multiple Ethernet interfaces as a result of a Linux kernel multicast socket routing bug.) UDP mode is activated when passing a single listening IP address and the builder object *second*:

.. code-block:: python

	collector = dfmux.DfMuxCollector("192.168.1.4", builder, [136])

If using multicast UDP for data transport, note that the mux system can deliver large numbers of UDP packets rapidly. If you see warnings about missed samples, you may need to increase the maximum size of the kernel UDP receive queue. On Linux, this can be accomplished by changing the value in ``/proc/sys/net/core/rmem_max``. On FreeBSD and Mac OS X, the maximum is in the sysctl ``kern.ipc.maxsockbuf``. A value of 5000000000 seems to work well.

On some versions of Linux with 128x DfMux firmware and multicast UDP for data transport, you will need to disable strict reverse-path validation in the kernel to take data. This can be accomplished by setting the sysctl ``net.ipv4.conf.all.rp_filter`` to 0. Depending on our system configuration, you may also need to set the corresponding per-interface sysctl (replace ``all`` with an interface name) corresponding to the network interface to which the DfMux boards are connected.

Lower data-loss rates with UDP can also often be achieved by setting the Qualityof-Service rules ("QoS") on your ethernet switch to respect DSCP indications (just look for the acronym).

Legacy Boards
_____________

This code can also be used to collect data from legacy boards with DAN firmware if you are so inclined by using the :py:class:`~spt3g.dfmux.LegacyDfMuxCollector` class in place of :py:class:`~spt3g.dfmux.DfMuxCollector`.

Data Structures
_______________

Frames generated by DfMuxBuilder contain two keys: "EventHeader" and "DfMux".

"EventHeader" is a :py:class:`~spt3g.core.G3Time` object containing the IRIG time of the first sample in the frame. If all the boards are synchronized correctly, this will also be the timestamp attached to all DfMux board samples.

"DfMux" is an object of type :py:class:`~spt3g.dfmux.DfMuxMetaSample`. This is a dictionary that maps board serial number to a :py:class:`~spt3g.dfmux.DfMuxBoardSamples` object. This in turn is a dictionary that maps readout module number (0-7) to a :py:class:`~spt3g.dfmux.DfMuxSample` object. This contains the IRIG timestamp for the data in its ``Timestamp`` member as well as a 128-element array of all the bolometer data in ``Samples``, stored with I and Q interleaved (so element 0 is channel 1/I, 1 is channel 1/Q, 2 is channel 2/I, etc.).

As an example:

.. code-block:: python

	channel2q = frame['DfMux'][frame['DfMux'].keys()[0]][0][3]

This retrieves data from the first board in the array, module 1, channel 2, modulation Q.

Housekeeping Data
~~~~~~~~~~~~~~~~~

DfMux board housekeeping is collected by the :py:class:`~spt3g.dfmux.Housekeeping.HousekeepingConsumer` class. It will query all of the boards in the most recent wiring map (see `The Wiring Map`) when a Housekeeping frame appears in the datastream, placing the results in the key ``DfMuxHousekeeping``. 

Housekeeping frames at fixed intervals can be generated using :py:class:`~spt3g.dfmux.Housekeeping.PeriodicHousekeepingCollector`. Note that collecting housekeeping information generates noise in detector timestreams and should be done only at times that you do not care about the data.

.. note:: 

	Housekeeping collecting can take up to a few seconds. If you are worried about pipeline stalls, you may want to run the housekeeping consumer in a subprocess (see ``G3Pipeline.Add()``).

The resulting data are stored in a :py:class:`~spt3g.dfmux.DfMuxHousekeepingMap` map, indexed by board serial number. This can be cross-correlated to the wiring map data. Mezzanines, modules, and channels stored in the elements are 1-indexed, matching the convention from pydfmux.

For ease of cross-correlation, there is a function :py:func:`~spt3g.dfmux.Housekeeping.HousekeepingForBolo` that can will look up the housekeeping information for a particular named bolometer.

.. code-block:: python

	hk = dfmux.HousekeepingForBolo(self.hkmap, self.wiringmap, 'Bolometer')

By default, this only returns information for the channel (notably containing the carrier amplitude and frequency). If you want the board, mezzanine, module, and channel information returned as a tuple, in that order, pass the keyword argument ``all_hk=True``.

Building Timestreams
====================

All analysis tools use data in the form of G3Timestreams, indexed by bolometer ID. Timestreams are typically stored in a Scan (see :doc:`frames`) frame, which is constructed from a wiring map and Timepoint frames using DfMuxCollator_.

The Wiring Map
~~~~~~~~~~~~~~

The wiring map, stored in a Wiring frame at the beginning of data taking, stores the mapping between bolometer ID and (Board Slot/Address, SQUID, Readout channel) tuples -- the information required to connect a :py:class:`~spt3g.dfmux.DfMuxMetaSample` object to bolometer samples. The wiring map is stored as the key ``WiringMap`` in an object of type :py:class:`~spt3g.dfmux.DfMuxWiringMap` in a Wiring frame. In almost all cases, this is inserted into the data stream by the :py:class:`~spt3g.dfmux.HardwareMapTools.PyDfMuxWiringMapInjector` module. This module is typically inserted as the first module following the :py:class:`~spt3g.dfmux.DfMuxBuilder` and takes a pydfmux hardware map as input (note: *not* a pydfmux session):

.. code-block:: python

	pipe.Add(dfmux.PyDfMuxWiringMapInjector, pydfmux_hwm=hwm)

DfMuxCollator
~~~~~~~~~~~~~

The :py:class:`~spt3g.dfmux.DfMuxCollator` class builds Scan frames (and timestreams) from Timepoint frames using the wiring map. Scan boundaries are signalled by the insertion of empty Scan frames into the data stream. When the :py:class:`~spt3g.dfmux.DfMuxCollator` object encounters a Scan frame, it will do the following:

	1) Accumulate all subsequent DfMux samples into two timestream maps, indexed by the bolometer IDs stored in the wiring map: ``RawTimestreams_I`` and ``RawTimestreams_Q``. Any samples for detectors not listed in the wiring map will be discarded. Accumulation ends with the next scan frame or the end of data processing, whichever comes first.
	2) Accumulate all scalar floating point numbers in the timepoint frames into timestreams with the same names. This is useful to store non-bolometer data such as telescope pointing.
	3) By default, FLAC compression is enabled for all bolometer timestreams and the source timepoint frames are discarded. These can be changed using the two arguments to the constructor of :py:class:`~spt3g.dfmux.DfMuxCollator`.

Empty scan frames can be inserted using a short Python module at appropriate boundaries. A trivial example is the :py:class:`~spt3g.dfmux.ScanTools.FixedLengthScans` module, which makes "scans" of some integer number of mux samples (by default, 1000 frames). In practice, you would want to break scans by GCP commands or analysis of telescope pointing.

.. code-block:: python

	pipe.Add(dfmux.PyDfMuxWiringMapInjector, pydfmux_hwm=hwm)
	pipe.Add(dfmux.FixedLengthScans, N=1000)
	pipe.Add(dfmux.DfMuxCollator)

Collecting data to a NetCDF file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :py:class:`~spt3g.dfmux.NetCDFDump` module takes timepoint frames and writes them to a NetCDF 3 file that can be opened using a variety of software packages, most notably KST, which will also monitor the file for updates. All sample points present in the wiring map are written to the output file with both I and Q demodulations, denoted by an ``_I`` or ``_Q`` suffix appended to the bolometer ID in the wiring map. In addition, a field called ``Time`` will be added containing the time of the sample (from the ``EventHeader`` key) in seconds since the UNIX epoch (Jan. 1, 1970). This time can be decoded using the python ``time`` module.

An example follows, including the addition of the wiring map from pydfmux and construction of the event builder:

.. code-block:: python

	pipe = core.G3Pipeline()
	builder = dfmux.DfMuxBuilder(len(hwm.query(pydfmux.core.dfmux.IceBoard).all()))
	collector = dfmux.DfMuxCollector("192.168.1.4", builder)
	pipe.Add(builder)

	# Insert current hardware map into data stream. This is critical to get the
	# channel -> board/module mapping needed to do anything useful with the data
	pipe.Add(dfmux.PyDfMuxHardwareMapInjector, pydfmux_hwm=hwm)

	pipe.Add(dfmux.NetCDFDump, filename=sys.argv[1])

This is contained in runnable form in ``dfmux/bin/ledgerman.py``.

Note that the version of KST installed from the default package repository under Ubuntu may not have support for reading NetCDF files produced by ledgerman. The version available from the KST PPA repository is compiled with NetCDF support (http://launchpad.net/~kst-plot/+archive/ubuntu/ppa).

Core Tools
==========

G3EventBuilder
~~~~~~~~~~~~~~

The :py:class:`~spt3g.core.G3EventBuilder` class implements an asynchronous frame builder based on frame objects delivered to its non-blocking ``AsyncDatum()`` call. When these arrive, the object calls the pure virtual method ``ProcessNewData()`` from a main thread. This method is responsible for assembling the data and eventually passing a complete frame to ``FrameOut()``, which will begin processing it in the pipeline. This is a C++-only abstract base class and is useful only when building a new data acquisition system.

G3TriggeredBuilder
~~~~~~~~~~~~~~~~~~

This C++-only class is the analog of G3EventBuilder for non-self-triggering systems (i.e. systems that poll for new data rather than streaming it). This can be used for once-every-N DAQ tasks like collecting housekeeping data.

ledgerman
=========

An example tool called ``ledgerman`` is included that collects data from the mux boards and writes it to a NetCDF file that can be read with kst. It is installed under ``bin`` in your build directory and will be available in your PATH if you have run ``env-shell.sh``.

.. code-block:: sh

	$ ledgerman /path/to/a/pydfmux/hardware/map.yaml output.nc

To see the frames as they go by:

.. code-block:: sh

	$ ledgerman -v /path/to/a/pydfmux/hardware/map.yaml output.nc

Like the other modules, you may see a few warnings about missing data immediately after it starts in the event that it starts collecting data midway through a sample. There should not be any warning messages after that.

