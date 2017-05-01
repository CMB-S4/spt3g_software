-----------
Quick Start
-----------

**Don't Panic**

.. contents:: Contents

How to Install
--------------

This depends on Boost and cmake, as well as the usual Python packages. Some additional packages (NetCDF, in particular) will activate optional components of the code if installed. You also need a C++11 compiler. This software is designed to run and work on a variety of operating systems (all Linuxes, Mac OS X, and FreeBSD) and architectures (at least 64-bit x86 and POWER).

Minimum versions:

	- GCC >= 4.7 or clang >= 3.3
	- Boost >= 1.48
	- cmake >= 2.6

Installing Dependencies on a Personal System
============================================

On Ubuntu/Debian, you can install the non-Python dependencies, including the optional ones, by doing:

.. code-block:: sh

	apt-get install cmake libboost-all-dev libflac-dev libnetcdf-dev libfftw3-dev libgsl0-dev

On RHEL-type systems (SL, CentOS, etc.), do this:

.. code-block:: sh

	yum install cmake netcdf-devel boost-devel flac-devel fftw-devel gsl-devel 

Note that on RHEL/SL versions before 7, you will need a newer compiler than ships with the OS. Please see the clustertools repository for a script in this case.

Installation on the Open Science Grid
=====================================

On an OSG or other system with OASIS configured, run this before anything else:

.. code-block:: sh

	eval `/cvmfs/spt.opensciencegrid.org/py2-v1/setup.sh`

This sets up a software environment with all the packages installed by yum, etc. above that you need for the SPT3G software environment, as well as a variety of standard cosmology and astrophysics tasks. You will obtain best results if you place the line above in your ``.bash_profile``. Do *not* put it in ``.bashrc`` and make *sure* that this is the *only* software installation set up in your bash profile. In particular, take care that there are no references to other python installations (Anaconda, etc.).

Compilation
===========

Having installed the appropriate dependencies, return to your checkout and run the following to build the software:

.. code-block:: sh

	mkdir build
	cd build
	cmake ..
	make

Passing ``-jN`` to ``make``, where N is the number of cores you wish to use during building, will speed up the process.

Once that is complete, you can use the ``env-shell.sh`` script in the build directory to set up the appropriate environment variables (PYTHONPATH, etc.):

.. code-block:: sh

	./env-shell.sh

Overview
--------

The large volume of SPT3G data, even for single observations, has forced some changes in the time-ordered-data processing workflow from previous processing to ensure that a minimum amount of data is in memory and being processed at any given moment. Typically, this minimum quantum of data is a left-right (or right-left) scan, which corresponds to the standard chunk size used in almost all filtering operations. You can of course also write code that runs on longer chunks of data, though this should be avoided where possible to avoid using too much memory. A short overview of the moving parts of the system appears below.

There are three main ingredients to data processing: frames, modules, and pipelines. Details on these topics can be found elsewhere in the manual, in particular in the :doc:`modules` and :doc:`frames` chapters; a brief overview is given here.


Frames
======

Frames (G3Frames) are generic data containers that behave like a python dictionary. They map arbitrary strings to arbitrary data. Here is an example:

.. code-block:: none

  In [31]: print frame
  Frame (Scan) [
  "ACUStatus" (spt3g.gcp.ACUStatusVector) => 3 elements
  "DfMuxHousekeeping" (spt3g.dfmux.DfMuxHousekeepingMap) => 37 elements
  "SourceName" (spt3g.core.G3String) => "RCW38"
  "GCPFeatureBits" (spt3g.core.G3VectorString) => 1 elements
  "RawBoresightAz" (spt3g.core.G3Timestream) => 386 samples at 190.783 Hz
  "RawBoresightEl" (spt3g.core.G3Timestream) => 386 samples at 190.783 Hz
  "RawTimestreams_I" (spt3g.core.G3TimestreamMap) => Timestreams from 1729 detectors
  "RawTimestreams_Q" (spt3g.core.G3TimestreamMap) => Timestreams from 1729 detectors
  "TrackerStatus" (spt3g.gcp.TrackerStatus) => 300 tracker samples from 21-Apr-2015:01:50:19.010000000 to 21-Apr-2015:01:50:22.000000000
  "Turnaround" (spt3g.core.G3Bool) => True
  ]

This frame contains information from a scan over RCW38 that you can access by the names in the first column, with a summary of their contents on the right. The (Scan) at the top is a description of the kind of data in the frame (e.g. Housekeeping data, a Map, a Scan, etc.)

The types of data you can store in the frame are containers that subclass G3FrameObject. These are listed in the manual for each Python module under the "Frame Objects" heading.

Modules
=======

A module is a python callable that does data processing. Modules are passed a frame and can inspect and modify it at will before the frame is passed along to the next module. An example of a module is doing poly filtering on a timestream. As an example of a very simple module:

.. code-block:: python

    def simplemod(frame):
        print(frame)

This prints the contents of the frame and does not modify it. As a more complex example, this would print the time at which a DfMux sample was recorded:

.. code-block:: python

    def printmuxtime(frame):
        print(frame['EventHeader'])

Modifying the frame also works like a dictionary. The following adds the number 5 to every frame:

.. code-block:: python

    def five(frame):
        frame['Five'] = 5

Much more detail is contained in the :doc:`modules` chapter of the documentation.

Pipelines
=========

A pipeline (G3Pipeline) is a sequence of modules. When the pipeline's Run method is invoked, it will run all modules in sequence for each frame in the data stream. Conceptually, it's nearly the same as a for loop. For example,

.. code-block:: python

    p = core.G3Pipeline()
    p.Add(dostuff)
    p.Add(dootherstuff)
    p.Run()

is equivalent to:

.. code-block:: python

    for frame in frames:
        dostuff(frame)
        dootherstuff(frame)

IO
==

Frames can be pickled and unpickled very quickly (1400 MB/s). Two special modules are provided (G3Reader and G3Writer) whose functions are to read and write frames to disk. This provides a full intermediate data format that can dump and restore the state of a pipeline to disk at any point. Something else equivalent to the above example is:

.. code-block:: python

    p = core.G3Pipeline()
    p.Add(dostuff)
    p.Add(core.G3Writer, filename='dump.g3')
    p.Run()

    p = core.G3Pipeline()
    p.Add(core.G3Reader, filename='dump.g3')
    p.Add(dootherstuff)
    p.Run()

You can also read files directly:

.. code-block:: python

    for frame in core.G3File('dump.g3'):
        dostuff(frame)

If for exploration you would like to load a file into memory the following idiom works.  Do not write code that relies on loading an entire file into memory or everything we've done was for naught.  This is just for poking at data:

.. code-block:: python

    frames = [fr for fr in core.G3File('thefilename.g3')]


Frame Objects
=============

Frames can store only objects that are subclasses of G3FrameObject or are plain-old-data (numerical scalars, booleans, strings). Notably, you cannot directly store python lists, tuples, or numpy arrays; container classes for these are provided, however. The primary driver for this is that the containers can be shared by C++ and Python code, which allows us to limit the amount of C++ to the cores of algorithms and preserve APIs between the two languages. This makes it much easier to write modules in C++ and Python interchangeably since both languages can access all the data products in the frame using the same interfaces.

The software provides both generic container classes (along the lines of a plain numpy array) and application-specific classes (such as ``G3Timestream``) that also contain metadata (for example, start and stop times and units). In general, code should use one of the purpose-specific objects, which makes sure that stored information has all the appropriate metadata attached.

Some classes that hold multiple instances of other datatypes have names starting with either G3Vector, which denotes a list/array, or G3Map, denoting a dictionary from strings to the named type. These names follow the C++ convention.

Classes containing large quantities of numbers (G3Timestream, G3SkyMap, G3VectorDouble) store their data contiguously in memory and implement the Python buffer protocol, which makes numpy operations on these classes behave with the same speed and semantics as on numpy arrays.

Experimental data is stored in one of the following application-specific clASSES:

* *G3Timestream* acts like a G3VectorDouble with attached sample rate, start time, stop time and units.
* *G3SkyMap* is a base class for actual maps of the sky, and includes units and projection information.
* *BolometerProperties* Stores the physical bolometer information like polarization angle and pointing offset.
* *DfMuxChannelMapping* Is used to map the string identifying a bolometer to its board/module/channel in the dfmux system.

A few notable generic containers when the standard ones are not appropriate:

* *G3VectorDouble* is a vector of doubles.  It acts like a numpy array of doubles.
* *G3MapString* acts like a dictionary that maps strings to strings
* *G3MapVectorDouble* acts like a dictionary that maps a string to a vector of doubles

Frame objects must be defined in both C++ and Python, which can be a bit daunting if you aren't familiar with C++.  If you *need* to add an extra member to a G3FrameObject subclass or need a new class, ask on the Slack channel and someone familiar with the C++ side of the software can help with it.

Units
=====

This software includes a units system that is meant to end wondering whether a given function takes radians or degrees as an argument, or whether a stored time is in milliseconds or seconds. The support code is accessible to both C++ and Python as part of the ``G3Units`` namespace (``core.G3Units.X`` in Python and ``G3Units::X`` in C++).

You should read the documentation on the :doc:`units` system.

Debugging Code
==============

Because of the step-by-step frame handling and callback system, debugging code requires a few more steps than usual.

To break into a debugger session at a certain point in the pipeline, you can use the ``spt3g.core.InjectDebug`` module.

Another common idiom is to insert a pipeline module that grabs data as it goes by for later examination, which lets you debug as though there were not callbacks. For example,

.. code-block:: python

    stuff = []
    def grabstuff(fr):
        if 'MyData' in fr:
            stuff.append(fr['MyData'])
    pipe.Add(grabstuff)

You can run the unit tests by running ``make test`` in the build directory, which is also a useful, though not sufficient, test that everything is working correctly -- expanding test coverage is always a praiseworthy activity.
