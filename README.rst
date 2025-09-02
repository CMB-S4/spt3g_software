.. image:: https://badge.fury.io/py/spt3g.svg
    :target: https://badge.fury.io/py/spt3g

.. image:: https://github.com/CMB-S4/spt3g_software/actions/workflows/cmake.yml/badge.svg
    :target: https://github.com/CMB-S4/spt3g_software/actions/workflows/cmake.yml

.. image:: https://github.com/CMB-S4/spt3g_software/actions/workflows/wheels.yml/badge.svg
    :target: https://github.com/CMB-S4/spt3g_software/actions/workflows/wheels.yml

About
-----

This repository contains the public portions of the analysis software developed by the SPT-3G collaboration. It includes the IO libraries, data acquisition code for McGill's DfMux readout boards, pipeline glue code, a binned map-maker, and a sky mock-observer, as well as defining all the APIs used by the project. For more detail on the intended processing architecture, please see the Quick Start chapter of the `documentation`_.

Except where otherwise noted, all files are distributed under the 2-clause BSD license. Acknowledgments and patches are appreciated.

.. _documentation: https://cmb-s4.github.io/spt3g_software/

Documentation
-------------

The main documentation for the software is in the docs folder. After building the software, you can build a pretty, searchable copy by running ``make docs``.

Dependencies
------------

This depends on pybind11 and cmake, as well as the usual Python packages. Some additional packages (NetCDF, in particular) will activate optional components of the code if installed. You also need a C++11 compiler. This software is designed to run and work on a variety of operating systems (all Linuxes, Mac OS X, and FreeBSD) and architectures (at least 64-bit x86 and POWER).

Minimum versions:

- GCC >= 5.0 or clang >= 3.4
- pybind11 >= 2.13
- cmake >= 3.12
- Python >= 3.7 (although pre-Python-3.8 support is best-effort)

On Ubuntu/Debian, you can install the non-Python dependencies, including the optional ones, by doing:

.. code-block:: shell

	apt-get install cmake libz-dev libbz2-dev liblzma-dev libflac-dev libnetcdf-dev

On RHEL-type systems (SL, CentOS, etc.), do this:

.. code-block:: shell

	yum install cmake netcdf-devel zlib-devel bz2-devel xz-devel flac-devel
	
If your system defaults to Python 2, but you wish to use Python 3, please do the following:

1. Install Python 3 *from the system package manager*
2. When you run cmake below, pass ``-DPython_EXECUTABLE=`which python3```

On any system, this software requires numpy and scipy (hard requirements), plus astropy and healpy (optional).

Setup on RHEL Systems with OASIS
--------------------------------

If using RHEL8 or RHEL9 and you have OASIS set up on your local system, run this before anything else:

.. code-block:: shell

	eval `/cvmfs/spt.opensciencegrid.org/py3-v5/setup.sh`

This will add a relatively new compiler and all necessary development packages (FLAC, GSL, etc) for compiling this software.


How to Build
------------

To build:

.. code-block:: shell

	cd spt3g_software
	mkdir build
	cd build
	cmake ..
	make

This will collect all of the necessary python tools into the ``build/spt3g`` directory, which can then be imported in python.  To set the appropriate python environment *without* installing the python package, use the shell script in ``build/env_shell.sh`` to run commands that know where to find the ``spt3g`` package and libraries:

.. code-block:: shell

	./env-shell.sh python my_script.py  # to run a python script
	./env-shell.sh ipython  # to start an ipython session

Alternatively, for users that only use a single build environment, set the following environment variables (e.g. in your ``.bash_profile`` file):

.. code-block:: shell

	export SPT3G_SOFTWARE_BUILD_PATH=path/to/spt3g_software/build
	export PYTHONPATH=$SPT3G_SOFTWARE_BUILD_PATH:$PYTHONPATH
	export LD_LIBRARY_PATH=$SPT3G_SOFTWARE_BUILD_PATH/lib:$LD_LIBRARY_PATH
	export PATH=$SPT3G_SOFTWARE_BUILD_PATH/bin:$PATH

Building the Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

You may find that you are missing some of the required packages for building the documentation.  To fix this, run the following commands:

.. code-block:: shell

	cd spt3g_software
	pip install -r doc/requirements.txt

To build the documentation in the build directory type:

.. code-block:: shell

	make docs

This will construct an html version of the documentation.  This builds the documentation in the build/docs folder.  Open build/docs/index.html in your favorite web browser.  You should at least read the quick start portion of the documentation before getting started.

Installation
------------

For various reasons it may be useful to install the software after building, instead of continuing to use it out of the build directory. Several CMake variables control how the software is installed:

* ``WITH_GZIP``, which defaults to ``TRUE``, is used to control whether the core library is built with support for gzip compression of G3 files.  Use ``-DWITH_GZIP=FALSE`` when calling ``cmake`` to disable.
* ``WITH_BZIP2``, which defaults to ``TRUE``, is used to control whether the core library is built with support for bzip2 compression of G3 files.  Use ``-DWITH_BZIP2=FALSE`` when calling ``cmake`` to disable.
* ``WITH_LZMA``, which defaults to ``TRUE``, is used to control whether the core library is built with support for lzma compression of G3 files.  Use ``-DWITH_LZMA=FALSE`` when calling ``cmake`` to disable.
* ``CMAKE_INSTALL_PREFIX``, which defaults to ``/usr/local`` is used as the root directory for installing all non-python components (header files, cmake export scripts, etc.).  This variable is frequently useful when installing into a python virtual environment.
* ``CMAKE_BUILD_PARALLEL_LEVEL`` is an environment variable (*not* a cmake option) used to control how many parallel processes are used to compile the shared libraries.  This option provides the same behavior as running ``make`` with the ``-j`` flag (e.g. ``make -j4``).

An uninstall target is also provided, so running ``make uninstall`` from the build directory should remove all files created by a previous ``make install``. 

Installation with Pip
---------------------

Use ``pip`` to install the python package.  Ensure that you use the appropriate options as necessary for your installation, e.g. ``--user`` or ``--prefix``.

For pre-built wheels hosted on `PyPI`_, available for most Linux x86_64, macOS x86_64 and macOS arm64 platforms, simply install the package without any additional options:

.. code-block:: shell

	pip install spt3g

The hosted wheels will include the necessary libraries (flac, etc) bundled with the package.  Otherwise, ensure that the dependency libraries are installed as explained above, and processed to one of the following steps.

To install the package from the github repo, run ``pip`` as usual (this may take a while, so consider setting the ``CMAKE_BUILD_PARALLEL_LEVEL`` environment variable):

.. code-block:: shell

	cd spt3g_software
	CMAKE_BUILD_PARALLEL_LEVEL=4 pip install -v .

By default this will create a directory called ``build`` in the repo and run the ``cmake`` build from there.  The build directory location can be changed by setting the ``BUILD_DIR`` environment variable, but keep in mind that ``pip`` requires that the build directory must be a path inside the repo file tree.
For development builds, use the ``--editable`` option to assemble the python package from the appropriate compiled extensions and python directories:

.. code-block:: shell

	cd spt3g_software
	CMAKE_BUILD_PARALLEL_LEVEL=4 BUILD_DIR=build pip install -v --editable .

An editable build adds references to the python directories to your python path, so that edits to library python files are immediately reflected in a fresh python session.

To pass arguments to the cmake build system, use the ``CMAKE_ARGS`` environment variable with arguments separated by spaces.  For example:

.. code-block:: shell

	cd spt3g_software
	CMAKE_ARGS="-DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_MODULE_PATH=/usr/local/share/cmake" pip install -v --prefix=/usr/local .

To run the test suite on the compiled package, you must have ``cmake``, and in particular the ``ctest`` utility, available on your path.  You must also know the location of the build directory where the cmake build was assembled (e.g. the value of ``$BUILD_DIR`` above).

.. code-block:: shell

	ctest --test-dir path/to/spt3g_software/build --output-on-failure

.. _PyPI: https://pypi.org/p/spt3g


Release Version Tracking
------------------------

Use git tags to keep track of release versions.  Tags should be of the form "v0.1.2" for release with major version 0, minor version 1 and patch version 2.
If such a tag is defined, cmake will populate the following outputs:

* A ``cmake/Spt3gConfigVersion.cmake`` file that contains the version number to be checked when including the Spt3g libraries in another cmake project
* A ``spt3g/version.py`` file containing VCS parameters for access in python and stored in PipelineInfo frames
* Add a ``SPT3G_VERSION`` compiler definition for accessing the version string in C++ code

Use the ``git archive`` command or the Python ``build`` package to export the source tree to a standalone archive.

Version Control Hygiene
-----------------------

The following is a brief overview of how to use git in a way that your collaborators will appreciate.

To initially check out the repository:

.. code-block:: shell

	git clone https://user@github.com/CMB-S4/spt3g_software.git

To update your checkout (the --rebase is important, especially if you have local changes):

.. code-block:: shell

	git pull --rebase

To send your changes back:

.. code-block:: shell

	git diff files_to_commit <- Examine this
	git commit files_to_commit
	git push

