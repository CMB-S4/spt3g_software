-------
File IO
-------

This framework has a native file format using the Cereal serialization library that can store all frames and maintain an exact copy of the data flowing through a pipeline. Since modules communicate with each other only through frames, interposing a serialized step in a pipeline allows it to be resumed from disk at any point. The file format is streaming (i.e. there are no headers) and has strong data integrity protections, including a CRC32-C checksum on each stored frame. The format is also architecture and endian safe.

Use of these files is through four main classes: G3Reader_, G3Writer_, G3MultiFileWriter_, and G3File_.

G3Reader
========

The G3Reader class reads a file from disk and inserts it into a pipeline when placed as the pipeline's first module. Once the file, or files, have been read completely, it will terminate the pipeline.

.. code-block:: python

	pipe = core.G3Pipeline()
	pipe.Add(core.G3Reader, filename='/path/to/file.g3')
	pipe.Add(core.Dump)
	pipe.Run()

To read multiple files, replace the path with a python iterable of paths:

.. code-block:: python

	pipe = core.G3Pipeline()
	pipe.Add(core.G3Reader, filename=['/path/to/file1.g3', '/path/to/file2.g3'])
	pipe.Add(core.Dump)
	pipe.Run()

All of the frames in ``file2.g3`` will appear to pipeline modules in order following all of the frames in ``file1.g3``.

G3Reader also supports transparent gzip decompression. If a file path ends in ``.gz``, it will be decompressed in memory as it is read:

.. code-block:: python

	pipe = core.G3Pipeline()
	pipe.Add(core.G3Reader, filename='/path/to/file.g3.gz')
	pipe.Add(core.Dump)
	pipe.Run()

In addition to reading files, G3Reader can be used to receive streaming data over the network; this is described in the :doc:`networkstreaming` section of the manual.

G3Writer
========

G3Writer is the output counterpart of G3Reader. All frames it receives as input get written to disk, as well as copied to the module's copied without changes. An example follows:

.. code-block:: python

	pipe = core.G3Pipeline()
	pipe.Add(dosomethings)
	pipe.Add(core.G3Writer, filename='/path/to/file.g3')
	pipe.Run()

All frames emitted by ``dosomethings`` will be written to ``/path/to/file.g3``. If this file is read in a different scipt using G3Reader, modules following the G3Reader will operate exactly as though they followed ``dosomethings``.

Like G3Reader, G3Writer supports gzip and will write compressed output if its output file name ends in ``.gz``.

G3MultiFileWriter
=================

G3MultiFileWriter is like G3Writer but, instead of writing to only one file, it splits its output over multiple files. This is useful, for example, when doing data acquisition or TOD streams to prevent having a single TB output file. In order to make these files independently readable, the most recent instances of any metadata frames (i.e. any frame of type other than Timepoint or Scan) will be prepended to every output file in the order in which they were originally seen.

The constructor of G3MultiFileWriter takes three arguments: a base file name, a file size limit, and, optionally, a file division algorithm. A typical invocation looks like the following and fills a directory with files of at most 1 GB:

.. code-block:: python

	pipe = core.G3Pipeline()
	pipe.Add(dosomethings)
	pipe.Add(core.G3MultiFileWriter, filename='/path/to/file-%02u.g3', size_limit = 1024**3)
	pipe.Run()

The base file name can be either a string or, for more complex naming schemes, a Python callable. The string uses printf-style formatting to substitute in a file sequence number. In the above case, the output files would be named file-00.g3, file-01.g3, file-02.g3, etc. It is also possible to pass a Python callable that returns a file name for more complex cases:

.. code-block:: python

	pipe = core.G3Pipeline()
	pipe.Add(dosomethings)
	pipe.Add(core.G3MultiFileWriter, filename=lambda frame, seq: '/path/to/file-%s-%d.g3' % (frame['SourceName'], seq), size_limit = 1024**3)
	pipe.Run()

Arbitrarily complex file division strategies can be employed using the divide_on argument. Like the base file name, this can be either static data -- an iterable of frame types -- or a Python callable. If passed an iterable of frame types, G3MultiFileWriter will always start a new file, even if the size limit is not yet reached, when it gets a frame of a type in the list. For example, to split the data into files of at most one GB, each of which includes data from at most one observation (each observation header frame will start a new file):

.. code-block:: python

	pipe = core.G3Pipeline()
	pipe.Add(dosomethings)
	pipe.Add(core.G3MultiFileWriter, filename='/path/to/file-%02u.g3', size_limit = 1024**3, divide_on=[core.G3FrameType.Observation])
	pipe.Run()

For more complex cases, you can also pass a callable as divide_on that takes a frame and returns ``True`` if a new file should be started without reference to the size limit and ``False`` otherwise. The following is equivalent to the above example, but uses a Python lambda function in place of the list:

.. code-block:: python

	pipe = core.G3Pipeline()
	pipe.Add(dosomethings)
	pipe.Add(core.G3MultiFileWriter, filename='/path/to/file-%02u.g3', size_limit = 1024**3, divide_on=lambda frame: frame.type in [core.G3FrameType.Observation])
	pipe.Run()

G3File
======

Unlike G3Reader and G3Writer, G3File is not a pipeline module. Instead, it is a python iterable that can be used to read frames interactively:

.. code-block:: python

	for frame in core.G3File('/path/to/file.g3'):
		print(frame)

It supports the same arguments as G3Reader.

File Format
===========

The file format used by these tools is a concatenation of serialized G3Frames. Each frame is serialized by the Cereal library using its portable binary archive mechanism in a simple endian- and word-length-neutral binary format. Each class is serialized using code in the class definition with fields stored sequentially on disk. The serialized form of each class is also versioned to allow later versions of the software to read earlier versions of the serialized classes.

Frame Structure on Disk
~~~~~~~~~~~~~~~~~~~~~~~

G3Frames are stored on disk with the following layout:

1. A 32-bit version code (currently set to 1) to describe the frame version
2. A 32-bit count of the number of objects stored in the frame.
3. The 32-bit type code.
4. A list of frame objects stored as:
	A. A string with the object key
	B. The blob for the frame object
5. A CRC32-C checksum

Blobs
~~~~~

Frame objects are stored in the frame through the intermediate structure of a blob. A blob is a portable binary archive of the frame object, serialized by a pointer to the base class, written to a binary buffer.

When a frame is read from disk, the blobs are read and stored in memory but not immediately deserialized. Instead, deserialization happens only lazily when the object is accessed. This optimizes file IO performance for data in which only a subset of keys are regularly being used and allows files to be read in which some keys are classes in unloaded libraries.

When a frame object is accessed (e.g. with []), the blob is deserialized. The original blob may, however, be retained in memory as well. As frame objects are only returned via const pointers, they are immutable for the life of the frame and so the original blob can be written back to disk later without expending CPU time to reserialize the class. Blobs are not retained for very large objects (> 128 MB) to lower memory consumption.

Miscellany
~~~~~~~~~~

Due to the streaming file format, the frames in one file can be appended to the frames in another file simply by using ``cat`` on the command line:

.. code-block:: bash

	cat file1.g3 file2.g3 > file-combined.g3


