Build System
------------

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

	#include <G3Frame.h>
	#include <pybindings.h>
	#include <boost/python.hpp>

	BOOST_PYTHON_MODULE(newthing)
	{
		boost::python::import("spt3g.core"); // Import core python bindings

		G3ModuleRegistrator::CallRegistrarsFor("newthing");
	}

This is sufficient for most uses (with "newthing" replaced by the name of the project).

Adding a C++ executable
~~~~~~~~~~~~~~~~~~~~~~~

You can also add C++ executables. Usually, there is not much reason to do this since everything is designed to be interacted with by Python. A few projects contain small standalone executables nonetheless, typically as test programs.

.. code-block:: cmake

	add_executable(newthingexec MyNewThingExecutable.cxx)
	target_link_libraries(newthingexec core newthing)
	list(APPEND SPT3G_PROGRAMS newthingexec)
	set(SPT3G_PROGRAMS ${SPT3G_PROGRAMS} PARENT_SCOPE)

The ``target_link_libraries`` command works as in `Adding a C++ library`_ above. The first command produces an executable named ``newthingexec`` that will be placed in the ``bin`` subdirectory of the build directory. The ``list`` and ``set`` commands inform other parts of the build system that this executable will exist, so that it can be included during installation. 

Mixing C++ and Python
=====================

If your project has both a C++ and a Python component, place the following into your ``__init__.py``:

.. code-block:: python

	from spt3g.core.load_pybindings import load_pybindings
	load_pybindings(__name__, __path__)

This (with no modifications) will merge the C++ and Python parts of the module into a single Python namespace.

