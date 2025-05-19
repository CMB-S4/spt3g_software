-------
Modules
-------

A module is a processing step for data. Modules operate on frames, and pass frames from their input to their output queues, applying some processing in between. A typical module will operate on one frame at a time, reading data (e.g. bolometer timestreams) from a single input frame, adding new data (e.g. a filtered version of the bolometer timestreams), and then pushing the frame to its output queue for the next module in the chain.

More complex manipulations are possible: the module can create new frames, throw frames away, operate on more than one frame at a time, or some combination of all three. Examples of how to do these things, as well as rationales, appear below.

.. contents:: Contents

Writing a module in Python
==========================

Modules can be written in Python in any one of three styles: as python functions, as subclasses of the base class core.G3Module, or as generic Python callables.

.. _function-modules:

Python Modules as Functions
___________________________

The simplest, and most common, case of a Python module is one that receives one frame, (potentially) modifies it, then passes it to the next module. This can be easily written as a simple Python function:

.. code-block:: python

	def simplemod(frame):
		print(frame)
	pipeline.Add(simplemod)

This prints its input frame to the console, then implicitly passes it on to the next module in the chain.

As frames behave like Python dictionaries, the same approach can be used to modify frames:

.. code-block:: python

	def five(frame):
		frame['Five'] = 5
	pipeline.Add(five)

The next module in the chain will now see a key named "Five" containing the number 5.

A full example of this doing something actually useful is to implement scan-by-scan poly-1 filtering of timestreams:

.. code-block:: python

	import scipy.signal

	def poly1(frame):
		if frame.type != core.G3FrameType.Scan:
			return
		outts = core.G3TimestreamMap()
		for i in frame['CalTimestreams']:
			outts[i.key()] = core.G3Timestream(scipy.signal.detrend(i.data(), units=i.data().units))
		outts.start = frame['CalTimestreams'].start
		outts.stop = frame['CalTimestreams'].stop
		frame['Poly1FilteredTimestreams'] = outts
	pipeline.Add(poly1)

This ignores non-scan frames and then creates a new timestream map containing a detrended version of the original, preserving the start and stop times and units. Right now, all the parameters of this processing step are hardcoded. It can be made configurable by the addition of keyword arguments:

.. code-block:: python

	import scipy.signal

	def poly1(frame, input='CalTimestreams', output='Poly1FilteredTimestreams'):
		if frame.type != core.G3FrameType.Scan:
			return
		outts = core.G3TimestreamMap()
		for i in frame[input]:
			outts[i.key()] = core.G3Timestream(scipy.signal.detrend(i.data()))
			outts[i.key()].units = i.data().units
		outts.start = frame[input].start
		outts.stop = frame[input].stop
		frame[output] = outts
	pipeline.Add(poly1, input='SomeOtherTimeStreams', output='OtherFilteredTimeStreams')

Note that the module does not modify the original timestreams in place. This is deliberate (and, in fact, modules in C++ are not even able to do this). The rationale here is that overwriting data in place:

  1. Makes it more confusing to trace the processing flow.
  2. Prevents some optimizations with file IO.
  3. Can create some causality paradoxes for certain kinds of data cached by modules.

Return values from Python modules
_________________________________

The examples above return ``None`` and so implicitly pass their input frame to the next module in the chain. Modules that need more control over data processing convey this by their return values:

  ``None``
    Passes input frame to the next module

  A G3Frame object
    Passes the return value to the next module **instead of** the input frame. This is usually used for the first module in a chain. The very first module has no data to work with and is responsible for generating it (see `The first module`_).

  An iterable of G3Frames
    Will insert the entire iterable (e.g. a Python list ``[]``) of frames into the input queue for the next module. This can be used to inject new data mid-processing, for example to read in calibration data, by returning a list containing both the input frame and a new one. Note that returning an empty list (``[]``) will cause the input frame to be dropped, which can be used to implement cuts. If the first module in the chain returns an empty list (``[]``), data processing will stop.

  Something with truth value (e.g. ``True`` or ``False``)
    A return value of ``True`` will cause the input frame to be passed to the next module and is equivalent to returning ``None``. Returning ``False`` will cause the input frame to be dropped and is equivalent to returning ``[]``. This can be used to implement cuts by returning the value of a conditional expression.
    
(Im)mutability of frame objects
_______________________________

Objects stored in a G3Frame are *immutable* -- once added to a frame, they should not be modified. If you want to change the contents of a frame object, delete it from the frame, modify a copy, and add the copy to be frame. (In C++, the compiler will prevent in-place modifications since ``G3Frame::Get()`` returns a const pointer, but Python doesn't have a concept of ``const`` so you just have to pay attention to the rules yourself and don't get the compiler's help). Copying the most common objects (e.g. G3TimestreamMaps) is designed to be particularly lightweight to facilitate this. 

The immutability allows three particularly useful features of the software:

 * Any pipeline module can cache an object (e.g. calibration data) by reference, without using excess memory. If later modules had the ability to modify frame objects, this could change the past in some sense by changing previous pipeline stages' view of the data after it was cached.
 * Serialization and deserialization of frame objects happens when frames are read from or written to disk and every time they pass through an MPI communicator or between members of a multiprocess group. To amortize this cost, G3Frame saves the original serialized copy of small (< 128 MB) frame objects when deserialized and restores that copy to disk/MPI/etc. when the frame is saved again rather than pointlessly reserializing the unchanged object. If an object read from disk/MPI is modified in place and the frame is re-serialized, it will save the *original version* of the object, leading to a great deal of puzzling down the line.
 * The provenance of any object in the frame can be tracked at a high level since the only time things can change is when they are added to or deleted from a frame. This makes tracing where data came from in the pipeline *vastly* easier.
 
**CAUTION**: An astute reader will notice that there are a few patterns that, while illegal in the strictest sense of the rules (they change frame-objects in place and a C++ compiler would reject them), will not create causality paradoxes (e.g. changing an object that a function has just added to a frame). Please proceed with great caution when writing such code.

For software connoisseurs: this is similar to, and exists for largely the same reasons as, the "single-static-assignment" architecture popular in modern compilers.

The first module
________________

The first module added to a ``G3Pipeline`` object is special: unlike all others, it does not act on input frames, since these frames cannot have come from anywhere. Instead, it is responsible for generating them. The ``G3Reader`` module is an example of this: it generates frames by reading them from disk.

Unlike all other modules, the first module will be passed ``None`` instead of a frame. This module then inserts data into the processing queue by returning new frames (see `Return values from Python modules`_). Data processing will stop when it returns an empty list (``[]``).

.. _class-modules:

Callable Objects as Functions
_____________________________

In addition to Python functions, any Python callable (anything that implements the ``__call__`` method) can be used as a processing module. This can be useful for processing steps that need to maintain state, such as a map making module that needs to keep its in-progress map between scans. Semantics and return values are the same as for Python functions (see `Return values from Python modules`_).

.. code-block:: python

	class NumberOfCalls(object):
		def __init__(self, Output='NCalls'):
			self.out = Output
			self.ncalls = 0
		def __call__(self, frame):
			self.ncalls = self.ncalls + 1
			frame[self.out] = self.ncalls
	pipeline.Add(NumberOfCalls, Output='Calls')

An alternative would be to subclass the ``core.G3Module`` class, which is more equivalent to the C++ mechanism but makes no practical difference at all, except that it will be automatically documented (see `Autodocumentation of modules`_). The only other difference is that the ``__call__`` method is renamed ``Process`` in this case:

.. code-block:: python

	class NumberOfCalls(core.G3Module):
		def __init__(self, Output='NCalls'):
			super(NumberOfCalls, self).__init__()
			self.out = Output
			self.ncalls = 0
		def Process(self, frame):
			self.ncalls = self.ncalls + 1
			frame[self.out] = self.ncalls
	pipeline.Add(NumberOfCalls, Output='Calls')

Autodocumentation of modules
____________________________

Preceding your module with the ``@core.indexmod`` decorator will allow the ``spt3g-inspect`` tool to list it. This should be used for processing steps designed for public use *only* rather than one-off functions for internal use in larger blocks of code.

For example:

.. code-block:: python

	@core.indexmod
	def printframe(frame):
		'''Print frame to console'''
		print(frame)

will produce the following output of ``spt3g-inspect``:

..

	--- Processing module: spt3g.example.printframe ---
	Print frame to console

All subclasses of ``core.G3Module`` (both in Python and C++) are automatically treated as though they were marked with this decorator.

.. _cxx-modules:

Writing a module in C++
=======================

The process of writing a processing module in C++ is similar to the Python one. C++ modules use a slightly different interface than Python; in particular, they behave like the callable object interface where all methods return lists.

A C++ module must inherit from the ``G3Module`` class. Data processing happens through the ``Process`` method, which takes two arguments: an input frame and an output queue. Output frames are pushed onto the queue; the semantics of this output queue are identical to those for Python processing modules returning lists.

.. code-block:: c++

	#include <G3Frame.h>
	#include <G3Module.h>
	#include <G3Data.h>
	#include <pybindings.h>

	class Five : public G3Module {
	public:
		void Process(G3FramePtr frame, std::deque<G3FramePtr> &out) {
			frame->Put("Five", G3IntPtr(new G3Int(5)));
			out.push_pack(frame);
		}
	};

	EXPORT_G3MODULE("exampleproject", Five, init<>(), "Adds five");

This example creates a module called ``Five`` that, like the earlier Python example, adds a key named ``Five`` to every frame. It is a part of the library "exampleproject" and will be accessible from Python as ``exampleproject.Five``.

Interaction with Python occurs through the ``EXPORT_G3MODULE()`` macro. The first two arguments are the library name and class to export. The third gives the arguments to the constructor (none, in this case). The fourth is the docstring visible for the class in Python. An example configurable version of the class follows:

.. code-block:: c++

	#include <G3Frame.h>
	#include <G3Module.h>
	#include <G3Data.h>
	#include <string>

	class Five : public G3Module {
	public:
		Five(std::string output = "Five") : output_(output) {}
		void Process(G3FramePtr frame, std::deque<G3FramePtr> &out) {
			frame->Put(output_, G3IntPtr(new G3Int(5)));
			out.push_pack(frame);
		}
	private:
		std::string output_;
	};

	EXPORT_G3MODULE("exampleproject", Five, init<optional<std::string> >(args("output")), "Adds five");

Here, the ``init<>`` arguments are modified to reflect that the configuration parameter is a string, that it is optional (leaving out the ``optional<>`` will make it mandatory), and that it maps to a Python keyword argument named "output". If your constructor takes multiple arguments, enclose the entire init section in parentheses to avoid preprocessor errors.

.. _pipesegments:

Pipeline Segments
=================

The use of pipeline segments allows you to have a canned collection of modules that can be added to a pipeline as though it were a single module. An example would be a pipeline segment that performs standard timestream filtering, which may be made of many separate modules but where specifying them individually would be tedious and prone to error.

A pipeline segment is defined by a Python function that is marked by the ``@core.pipesegment`` decorator and takes a pipeline as its first argument. For example:

.. code-block:: python

	@core.pipesegment
	def standardfiltering(pipe, input='CalTimestreams', output='OutTimestreams'):
		'''
		This is the standard timestream filtering used for 2016 data
		'''

		pipe.Add(analysis.PolyFilter, input=input, order=1)
		pipe.Add(analysis.MaskedHighPassFilter, ell=3000)

	pipe.Add(standardfiltering, output='FilteredTimestreams')

By default, the ``core.pipesegment`` decorator will introspect these functions by running them against a fake pipeline object. This information about what the segment does is then automatically appended to the docstring for the segment. This makes it easy for a user to discover what your wrapper does in a way that cannot become inconsistent with documentation. If your pipeline segment has side effects (e.g. opening files) or cannot be run with its default arguments, you may wish to disable this behavior by passing the ``autodoc=False`` keyword argument to the decorator.

Advanced Techniques: Buffering Data
===================================

Modules that need to work on granularilty coarser than a scan (e.g. notch filtering) can **buffer** frames. This can be implemented using the Python callable interface. For example:

.. code-block:: python

	class Buffered(object):
		def __init__(self):
			self.buffer = []
		def __call__(self, frame):
			if len(self.buffer) < 5:
				# Add to buffer and move to the next scan
				self.buffer.append(frame)
				return []

			# Now we have 5 frames queued up
			dostuffwithfivescans(self.buffer)

			# Clear buffer and send these frames onward
			returnval = self.buffer
			self.buffer = []
			return returnval

This implements a processing step that works on five scans at a time. From the perspective of a module either before or after this one in the chain, nothing unusual happens: frames appear in order one at a time in both cases. When ``__call__`` returns an empty list, the pipeline goes back to the first module to get a new frame instead of continuing to the next. These accumulate inside the internal queue of ``Buffered`` until there are five scans present. At that point, they are processed as a group and then moved to the output queue. When the pipeline sees five frames in the output queue, it will call the next module five times, with each frame in sequence. Once that is complete, it will then go back to the first module for new frames.

Caching Previous Data
=====================

The previous example can be adapted to cache data from previous frames when that data is required to operate on current frames.  This is useful for caching e.g. calibration data.

.. code-block:: python

        class Caching(object):
                def __init__(self, calibrator_key='CalibratorResponse'):
		    self.calkey = calibrator_key
		    self.cal = None
		def __call__(self, frame):
		    if self.calkey in frame:
		        self.cal = frame[self.calkey]
			return
		    if frame.type == core.G3FrameType.Scan:
		        # Do something here

There are two convenience methods for caching data and passing it into a module.  This is useful for writing a function that requires (for example) a ``BolometerPropertiesMap``, but otherwise does not require a class.  This is best illustrated with an example

.. code-block:: python

    @core.cache_frame_data(type=core.G3FrameType.Scan, bolo_props='BolometerProperties')
    def FlagSomeStuff(frame, flag_key='Flags', bolo_props=None):
        pass

The decorator ``@core.cache_frame_data`` will ensure that the bolometer properties are passed to ``FlagSomeStuff`` in the ``bolo_props`` kwarg.  The (required) ``type`` argument specifies that ``FlagSomeStuff`` only runs on frames of type ``core.G3FrameType.Scan``.  When called with an abritrary ``kwarg='FrameKey'`` pair, ``core.cache_frame_data`` caches the most recent instance of ``'FrameKey'`` in any frame, and passes it to the decorated function under the keyword argument ``kwarg``.  ``core.cache_frame_data`` can cache multiple keys from multiple frames.  For example, if one wanted to cache the calibrator singal to noise ratio (stored as ``'CalibratorResponseSN'``) as well as the bolometer properties, and pass it to a function with the keyword argument `calsn`, one would wrap that function with the decorator ``@core.cache_frame_data(type=core.G3FrameType.Scan, calsn='CalibratorResponseSN', bolo_props='BolometerProperties')``.  This would be useful if one wanted to apply different calibrator SNR thresholds to different bands, for example.

``core.scan_func_cache_data`` is a special case of ``core.cache_frame_data`` where ``type`` is set to ``core.G3FrameType.Scan``.

Pipelines
=========

Modules are connected to one another by a pipeline object, of which there is currently one implementation: G3Pipeline. Any pipeline has two interesting methods, ``Add`` and ``Run``.

Pipeline.Add
____________

The ``Add()`` method adds a module to the pipeline immediately following the last added module. It accepts any of the module types described above, as well as pipeline segments. For classes (either C++ or Python), it can accept either an instance of the module or the class. If passed a class, keyword arguments following the class are passed to the class constructor. The following two pieces of code are equivalent:

.. code-block:: python

	pipe = G3Pipeline()
	pipe.Add(core.G3Reader(filename="test.g3"))

.. code-block:: python

	pipe = G3Pipeline()
	pipe.Add(core.G3Reader, filename="test.g3")

For pipeline segments, only the second syntax works. As a result, the second syntax is generally preferred, as it can be used uniformly for all objects that can be passed to ``Add()``. Additionally, only the second syntax will record configuration information (see G3PipelineInfo_).

``Add()`` accepts a special keyword argument (``name``) that can be used to set the name of a module or segment in the output of run profiling (see below). If unspecified, it defaults to the name of the class or function, with slashes indexing modules added by pipeline segments.

If the ``subprocess`` argument to ``Add()`` is set to True, the module passed will be run in a python subprocess using the multiprocessing framework. Note that this does *not* work yet for segments.

Pipeline.Run
____________

The ``Run()`` method runs the pipeline until completion (see `The first module`_). It takes one optional keyword argument (``profile``). If set to ``True``, it will print out the amount of system and user time spent in that module during processing after completion.

G3PipelineInfo
______________

G3Pipeline will automatically insert information about its configuration into the data stream by internally emitting a PipelineInfo frame containing a timestamped G3PipelineInfo object with the following information:

- Version control information (branch, revision number, source URL, version name if any, presence of local diffs, etc.) reflecting the software currently running.
- The user and host running the software.
- The configuration of all modules and/or segments added to the pipeline.

This information is added immediately following the first added module or segment. If the first frame in the data stream at this point is already a PipelineInfo frame (or a PipelineInfo frame occurs in one of the first few frames, with only metadata frames before it), the G3PipelineInfo object described above will be added to it; otherwise, a new PipelineInfo frame with the object is prepended to the data stream.

Within some limits imposed by Python (related to lambda functions, most notably), calling ``repr()`` on a G3PipelineInfo object (or a G3Pipeline object) will yield an executable Python script reflecting the exact modules and configuration used to produce the data. To within the mentioned limitations, this script can be rerun to exactly reproduce stored data; it can also be inspected to learn the configuration of the data's source pipeline[s] and thus the processing that produced it.

The G3PipelineInfo ``.Run()`` method provides a convenient way of rerunning the pipeline configuration within it, and the ``.modules`` attribute is a list of G3ModuleConfig objects with dict-like access to the arguments provided to each pipeline module.

Limitations:

- The content of functions defined inline in a script (either by ``def`` or ``lambda``), as opposed to functions defined in an imported Python module, will not appear in the output, though options will. Inline functions defined by ``def`` will at least give the name of the function.
- Options passed to pre-instantiated modules will not be stored. Only options passed in ``pipe.Add()`` will be recorded. For example, ``pipe.Add(core.G3Reader, filename="test.g3")`` will fully record its arguments, but ``pipe.Add(core.G3Reader(filename="test.g3"))`` will not. Prefer the syntax that records options unless you have a compelling reason to do something else.
- A G3Pipeline created in C++ will not record configuration; only G3Pipelines created in Python will.
- If the code used to run the pipeline initially has changed, rerunning the configuration from the G3PipelineInfo object will produce different results.  Use the version control information stored in the object to ensure you are running the correct version of the software.
