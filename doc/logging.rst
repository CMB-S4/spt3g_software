Logging
-------

This software includes a logging framework similar to the standard python logging that can be used from both Python and C++. There are seven log **levels** (below) that can be set for each logging **unit**. The default log level is ``log_notice``.

 log_trace
   Messages of the type "At line 75" that no one really ever wants to see except in the deepest throes of debugging. In C++ code, these messages are compiled out unless ``CMAKE_BUILD_TYPE`` is set to Debug.
 log_debug
   Messages of the type "Configuration parameters: Foo Bar" that almost no one really ever wants to see except while debugging. In C++ code, these messages are compiled out unless ``CMAKE_BUILD_TYPE`` is set to Debug.
 log_info
   Messages slightly more interesting than ``log_debug`` (for example: "Opening file X.h5"). Always available, but not shown by default.
 log_notice
   Messages that do not indicate a potential problem but that it is good for humans to see. For example, "Starting new output file X". Shown by default, so please limit the use of messages this level and above.
 log_warn
   Messages that indicate a potential problem such as "No data from boards in 30 seconds".
 log_error
   Serious recoverable data problem. For example, "Scan does not have input key". Printed in bold red so you know it is serious.
 log_fatal
   An unrecoverable problem such as a bad configuration. Throws an exception (RuntimeError) and stop all data processing.

Setting Log Levels
==================

The global logging level can be set to any of the above values as follows:

.. code-block:: python

	core.set_log_level(core.G3LogLevel.LOG_WARN)

The log level for an individual section of the code can be set as well, which is useful if you want to increase (or decrease) the verbosity of a specific module for debugging.

.. code-block:: python

	core.set_log_level_for_unit(core.G3LogLevel.LOG_WARN, 'NoisyThing')

Using Logging from Python
=========================

From Python, you can emit log messages using the ``core.log_*`` family of functions. For example:

.. code-block:: python

	core.log_error('Everything is on fire')

By default, all log messages from Python are assigned to the log unit "Python". If you pass a second argument, you can set the log unit to something else:

.. code-block:: python

	core.log_error('Everything is on fire', 'FireWarden')

Then you could selectively disable "FireWarden" if you don't want to see such messages by setting FireWarden's log level to log_fatal (see `Setting Log Levels`_).

Using Logging from C++
======================

The C++ interface to logging mirrors the Python one. For a C++ version of the FireWarden above, we can do the following:

.. code-block:: c++

	log_error("Everything is on fire");

In C++, the log functions know how to do ``printf()`` formatting. For example:

.. code-block:: c++

	log_error("%d things are on fire", 15);

Log units in C++ are set for a particular scope rather than in the call to the log function. This is done using the ``SET_LOGGER()`` macro. For example:

.. code-block:: c++

	void WardFire(void) {
		SET_LOGGER("FireWarden");

		log_error("%d things are on fire", 15);
	}

A common idiom for this is to put ``SET_LOGGER()`` into the class definition so that it applies to all log messages emitted by class member functions. As an example:

.. code-block:: c++

	class FireWarden : public G3Module {
	public:
		void Process(G3FramePtr frame, std::deque<G3FramePtr> &out) {
			log_error("%d things are on fire", frame->Get<G3Int>("ThingsOnFire")->value);
			out.push_back(frame);
		}
	private:
		SET_LOGGER("FireWarden");
	};

