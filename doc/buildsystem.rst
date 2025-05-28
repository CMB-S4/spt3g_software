Build System
------------

.. contents:: Contents
   :local:

CMake Overview
==============

We use cmake as a build system, which provides a number of nice features for finding dependencies and managing the build process. cmake generates Makefiles (or XCode projects) from files called ``CMakeLists.txt`` in each directory. cmake in general is designed to perform so-called "out-of-tree" builds, which keep the source directories unmodified by object files. The process for doing this is something like this:

.. code-block:: sh

	mkdir build
	git clone https://github.com/CMB-S4/spt3g_software src
	cd build
	cmake ../src
	make

A number of variables can be set on the command line when cmake is run that control the build. The syntax for these options is ``cmake -DVARIABLE=value srcdir``.

  CMAKE_BUILD_TYPE
    This can be set to either ``Release`` or ``Debug``. Setting it to ``Release`` will cause the compiler to optimize the code more, making it substantially faster at the expense of increased build times and removal of some debugging information.

  BUILD_PROJECTS
    This variable can be set to a semicolon-separated list of projects to build, to allow building only a subset of the projects which are present. For example, specifying ``-DBUILD_PROJECTS='gcp;dfmux'`` will result in the ``core``, ``gcp``, and ``dfmux`` projects being built. The ``core`` project must always be built, so if this variable is set to a list which does not contain it, it will be added. Other project dependencies are not automatically detected, so use this feature only if you are certain that you understand exactly which projects you need. Note that the list will usually need to be quoted to avoid the shell interpreting the first semicolon as the end of the command.

Adding a Project
================

To add another project to the repository (corresponding to a python module, e.g. ``from spt3g import newthing``), make a new directory at the root of the repository. This directory must contain a file called CMakeLists.txt.

Adding Python code
~~~~~~~~~~~~~~~~~~

To include Python code, include a line like the following in your CMakeLists.txt:

.. code-block:: cmake 

	execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR}/python ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/newthing)

Here, ``newthing`` is the name of your project (which **must** match the directory name) and ``python`` is the name of the directory inside your project containing the Python source code. If your project contains only Python code, this need not be a subdirectory and you can remove the "/python".

Adding a C++ library
~~~~~~~~~~~~~~~~~~~~

To add a C++ component to your project, add some lines like the following:

.. code-block:: cmake

	add_library(newthing SHARED
		src/MyNewThing.cxx src/python.cxx
	)
	target_link_libraries(newthing core ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

This builds a library called ``newthing.so`` from the two given source files and links it to the spt3g core library and the Boost and Python libraries (which are mandatory). Typically C++ source files are in a directory called ``src``. Header files that you want visible from other projects must be placed in a directory called ``include/newthing``.

Every C++ library must contain a file declaring the library to Python. This file is usually named ``python.cxx`` and has contents like the following:

.. code-block:: c++

	#include <pybindings.h>

	SPT3G_PYTHON_MODULE(newthing)
	{
		py::import("spt3g.core"); // Import core python bindings
	}

This is sufficient for most uses (with "newthing" replaced by the name of the project).

Adding a C++ executable
~~~~~~~~~~~~~~~~~~~~~~~

You can also add C++ executables. Usually, there is not much reason to do this since everything is designed to be interacted with by Python. A few projects contain small standalone executables nonetheless, typically as test programs.

.. code-block:: cmake

	add_spt3g_executable(newthingexec MyNewThingExecutable.cxx)
	target_link_libraries(newthingexec core newthing)
	list(APPEND SPT3G_PROGRAMS newthingexec)
	set(SPT3G_PROGRAMS ${SPT3G_PROGRAMS} PARENT_SCOPE)

The ``target_link_libraries`` command works as in `Adding a C++ library`_ above. The first command produces an executable named ``newthingexec`` that will be placed in the ``bin`` subdirectory of the build directory. The ``list`` and ``set`` commands inform other parts of the build system that this executable will exist, so that it can be included during installation. 

Adding tests
~~~~~~~~~~~~

Tests can be written in either Python or C++. Some tests must be written in one language  in order to test interfaces specific to it; otherwise, most tests are currently written in python. 

The simplest way to run the full set of tests is by executing ``make test``. This does not allow for much flexibility, however, so in cases where more control is desirable, one should run tests using the ``ctest`` driver tool directly. Commonly useful options are ``ctest --output-on-failure`` which will show a test's output when it fails, which is frequently useful for understanding what the failure was in order to fix it, and ``ctest -R <regex>`` which will run any tests any part of whose name is matches the given regular expression, which is handy for running just a particular test to debug it, without having to wait while the entire test suite runs each time. 

Python Tests
^^^^^^^^^^^^

Python tests should be placed in a ``tests`` subdirectory of the project. Each test must then also be declared in the project's ``CMakeLists.txt``, so that ``cmake`` will know to include it in the list of tests to be run by ``ctest`` or the ``test`` build target. This is done by using the ``add_spt3g_test`` macro:

.. code-block:: cmake

	add_spt3g_test(test_foo)

will add a test which is implemented in ``tests/test_foo.py``. 

The contents of a Python test script can be anything; the script is simply run, and if its exit value is 0, it is considered to have passed. Any non-zero exit status will be taken to indicate a failure. The simplest mechanism to do tests is to just write code with ``assert`` statements which check that properties of interest hold.

C++ Tests
^^^^^^^^^

C++ tests consist of one or more implementation files which declare tests, organized into test groups. The implementation files for a test are linked together into a test executable. 

Like Python tests, C++ tests must be declared in the project's CMake script, which is done using the ``add_spt3g_test_program`` command:

.. code-block:: cmake

	add_spt3g_test_program(test
	                       SOURCE_FILES
	                         ${PROJECT_SOURCE_DIR}/my_test.cpp
	                       USE_PROJECTS core)

The first argument is the name of the test executable, which will be prefixed with the project name. Several implementation files may be listed after ``SOURCE_FILES``, and the arguments after ``USE_PROJECTS`` indicate which projects the executable depends on, so suitable compiler options will be generated to give access to those projects' header paths and to link against their libraries. 
Arbitrary labels can also be associated with a test by passing them after the ``TEST_LABELS`` argument.

Typically, each implementation file defines one test group, but multiple implementation files may redeclare and contribute to the same test group. It is also possible to place multiple test groups in one translation unit by isolating each in its own namespace. Each test implementation file should include the ``G3Test.h`` header to have access to the test infrastructure definitions.

A test group is declared using the ``TEST_GROUP`` macro:

.. code-block:: c++

	TEST_GROUP(MyTests)

Individual tests are then defined using the ``TEST`` macro, followed by a function body which does the work of the test:

.. code-block:: c++

	TEST(Test1){
		Num::InitializeNumbers();
		auto n5 = Num::Get(5);
		auto n7 = Num::Get(7);
		ENSURE(n5 < n7, "5 should be less than 7");
	}

The argument to the ``TEST`` macro is the name of the test, which will then have a fully qualified name derived from its test group: ``MyTests::Test1``. 

Since multiple C++ tests can run in the same executable, it is poor form to use ``assert``, ``exit``, or some other mechanism which can stop the whole process before other tests can run. Tests indicate failure by throwing an exception, but for convenience and readability, particularly of failure messages, a set of macros are provided. The simpest is ``ENSURE``, which takes a predicate to be tested, and optionally a message to be shown if the predicate evaluates to false. An example is shown above, and if that test fails, the output produced would look similar to the following:

.. code-block:: none

	MyTests::Test1: /some/path/my_test.cpp:50: n5 < n7: 5 should be less than 7
	FAIL

The ``FAIL`` macro can be used when a test has reached a point in its control flow which indicates failure without any further condition needing to be checked. This is particularly useful for ensuring that exceptions are or aren't thrown in correct places:

.. code-block:: c++

	TEST(Exceptions){
		try{
			some_func();
		}
		catch(...){
			FAIL("some_func must not throw exceptions");
		}
		
		try{
			other_func(bad_val);
			FAIL("other_func must throw an exception when passed bad_val");
		}
		catch(...){}
	}

There is also the ``ENSURE_EQUAL`` macro, which specifically checks to expressions for equality, and produces a detailed error message if they are not:

.. code-block:: c++

	TEST(Equality){
		int a=4, b=5;
		ENSURE_EQUAL(a,b,"a and b should be the same");
	}

which outputs:

.. code-block:: none

	MyTests::Equality: my_test.cpp:19: ENSURE_EQUAL(a, b): 4 != 5: a and b should be the same

Mixing C++ and Python
=====================

If your project has both a C++ and a Python component, place the following into your ``newthing/python/__init__.py``:

.. code-block:: python

	from .._libnewthing import *

This will merge the C++ and Python parts of the module into the ``newthing`` Python namespace.

