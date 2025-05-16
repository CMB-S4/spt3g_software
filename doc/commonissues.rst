---------------
Common Problems
---------------

This contains a list of common problems that people have run into and solutions to them. If you run into others, and an FAQ entry here is a more appropriate place than text in the main part of the manual, please feel to add them to this page.

Setting Values of a FlatSkyMap
------------------------------

Because the implementations of ``__setitem__`` and ``__getitem__`` in G3SkyMap have to handle both 1D and 2D semantics), the FlatSkyMap class has different slicing semantics than a numpy array. To slice it as though it were a numpy array, you can cast it to one:

.. code-block:: python

	numpy.asarray(your_flat_sky_map)[:] = the_numpy_array_you_are_assigning

Because sky maps, like timestreams, implement the Python buffer protocol, there is no performance or memory penalty for these kinds of operations.


Error Message Decryption
------------------------

Argument Errors
~~~~~~~~~~~~~~~

If you are adding a module to a pipeline in python, and if the module is written in C++, when you misspell a keyword argument, you may get an unhelpful message like the following:

.. code-block:: none

	TypeError: some_complicated_function(): incompatible function arguments. The following argument types are supported:
	    1. ...

	Invoked with: ...

This message arises because Python and C++ have different semantics for overloaded functions, preventing pybind11 from succesfully disambiguating a misspelled or otherwise wrong argument from an attempt to execute a different function entirely. If you see errors like the above, check your arguements for typos in keywords, number of positional arguments, and types.

Misbehaving Compilers
---------------------

1) Clang 3.6.0 has some bugs with the std::unordered_map implementation that our code encounters.

