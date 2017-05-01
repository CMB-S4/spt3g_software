-----------------
Network Streaming
-----------------

The spt3g.core module provides facilities for consuming data from a network host and serving it, through the ``core.G3Reader`` and ``core.G3NetworkSender`` classes. Both exchange frames in the same wire format used by files read by ``core.G3Reader`` over TCP sockets.

While communication is unidirectional from the ``core.G3NetworkSender`` (or netcat) to a ``core.G3Reader`` (or netcat), either end of the connection can be a "server" or "client" in the sense of which end of the connection is expected to be listening for the other.

.. contents:: Contents

G3NetworkSender
_______________

G3NetworkSender_ can operate in two modes: as a reliable transfer mechanism that connects to a remote G3Reader or as a streaming data source to which remote clients can connect to get the data currently passing through the pipeline.

Connecting to a Remote Reader
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the simplest mode, a remote G3Reader (see below) is set up to listen for an incoming connection and G3NetworkSender_ is pointed at it as follows:

.. code-block:: python

    pipe.Add(core.G3NetworkSender, hostname='name.of.remote.host', port=4536)

This will send every frame passing through the pipeline to the remote host. Note that the output frames are transferred from a second thread, which builds its own unbounded buffer. In the event that file IO is much faster than network IO, this may have the effect of transiently reading the entire input into memory. (This can and should be fixed in a future version of the software)

Functioning as a Streaming Server
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the standard operational mode of G3NetworkSender_. When operating as a streaming server, remote clients connect to it and receive the most recent of any metadata frames (Calibration, Wiring, etc. -- anything that is not a Timepoint or Scan frame), followed by all passing data beginning at the time of connection. This is designed for use, for example, in monitoring data acquisition.

Streaming server mode is established by using the special hostname ``*``. To create a streaming server:

.. code-block:: python

    pipe.Add(core.G3NetworkSender, hostname='*', port=4536)

This can be connected to on the remote end by doing:

.. code-block:: python

    pipe.Add(core.G3Reader, filename='tcp://name.of.streamer.host:4536')

By default, this will send every frame starting at the instant of connection. For a monitoring system with potentially unreliable or slow data links, setting the ``max_queue_size`` option will limit the growth of buffers, dropping Timepoint or Scan frames in the event of network latency rather than buffering a number of frames larger than the ``max_queue_size`` parameter.

For example, the following will only allow the remote end to fall ten frames behind before it starts dropping data. When the remote end comes back from a potential transient disruption, it will then be in the present rather than slowly clearing a backup using memory on the server. Note that metadata frames (Housekeeping, Wiring, etc.) important for maintaining the state of the pipeline will *never* be dropped, no matter what the setting of ``max_queue_size`` is.

.. code-block:: python

    pipe.Add(core.G3NetworkSender, hostname='*', port=4536, max_queue_size=10)

Throttling output data
~~~~~~~~~~~~~~~~~~~~~~

An additional class called G3ThrottledNetworkSender wraps G3NetworkSender_ but with a data reduction step that sends only every Nth of some list of frame types. For example, the following sends all non-Timepoint frames and every 10th Timepoint frame with a configuration otherwise identical to the example in the previous section:

.. code-block:: python

    pipe.Add(core.G3ThrottledNetworkSender, hostname='*', port=4536, max_queue_size=10, frame_decimation = {G3FrameType.Timepoint: 10})

Other frame types can be added by appending them to the dictionary. All instances of frame types not occurring in the dictionary will be sent.

Using G3Reader Over the Network
_______________________________

G3Reader, in addition to opening files, can read from network sockets by passing it a URL of the form ``tcp://host:port`` instead of a file path.

To connect to the streaming server examples above, for example, you can add a G3Reader with the following configuration:

.. code-block:: python

    pipe.Add(core.G3Reader, filename='tcp://name.of.streamer.host:4536')

Real files and network sockets can be mixed. For example, if you want to stream real-time data from the DAQ system with some calibration information prefixed:

.. code-block:: python

    pipe.Add(core.G3Reader, filename=['/path/to/some/cal/data', 'tcp://name.of.streamer.host:4536'])

If the remote end closes the connection (as happens when the pipe in which G3NetworkSender_ is running finishes), G3Reader will interpret that as the end of the file and either move onto the next file or exit, as appropriate.

Like G3NetworkSender_, G3Reader also has a listen mode in which it will wait for a connection on a given port from an external source instead of connecting to a remote host itself. The semantics of its operation are otherwise identical to the connect-to-remote-host mode. Like G3NetworkSender_, this mode is triggered by the special host name ``*``. To make G3Reader wait for a remote connection on port 4536:

.. code-block:: python

    pipe.Add(core.G3Reader, filename='tcp://*:4536')


