About
-----

This repository contains the public portions of the analysis software developed by the SPT-3G collaboration. It includes the IO libraries, data acquisition code for McGill's DfMux readout boards, and pipeline glue code and defines all the APIs used by the project. For more detail on the intended processing architecture, please see the Quick Start chapter of the `documentation`_.

Except where otherwise noted, all files are distributed under the 2-clause BSD license. Acknowledgments and patches are appreciated.

.. _documentation: https://cmb-s4.github.io/spt3g_software/

Documentation
-------------

The main documentation for the software is in the docs folder. After building the software, you can build a pretty, searchable copy by running ``make docs``.

Dependencies
------------

This depends on Boost and cmake, as well as the usual Python packages. Some additional packages (NetCDF, in particular) will activate optional components of the code if installed. You also need a C++11 compiler. This software is designed to run and work on a variety of operating systems (all Linuxes, Mac OS X, and FreeBSD) and architectures (at least 64-bit x86 and POWER).

Minimum versions:

	- GCC >= 4.7 or clang >= 3.3
	- Boost >= 1.48
	- cmake >= 3.1
	- Python >= 2.7 (although pre-Python-3 support is best-effort)

On Ubuntu/Debian, you can install the non-Python dependencies, including the optional ones, by doing:

.. code-block:: shell

	apt-get install cmake libboost-all-dev libflac-dev libnetcdf-dev libfftw3-dev libgsl0-dev

On RHEL-type systems (SL, CentOS, etc.), do this:

.. code-block:: shell

	yum install cmake netcdf-devel boost-devel flac-devel fftw-devel gsl-devel
	
If your system defaults to Python 2, but you wish to use Python 3, please do the following:

	1. Install Python 3 *from the system package manager*
	2. Make sure the python-3 version of the Boost library is installed (on Ubuntu, this is part of the standard boost-python package referenced above)
	3. When you run cmake below, pass ``-DPYTHON_EXECUTABLE=`which python3```


Setup on RHEL6
--------------

Note that on any RHEL6 system, you will need a newer compiler than ships with the OS. Please follow whatever directions apply at your site to achieve this. Alternately, if you have OASIS set up on your local system, run this before anything else (also works on several other operating systems, including RHEL7):

.. code-block:: shell

	eval `/cvmfs/spt.opensciencegrid.org/py3-v2/setup.sh`


How to Build
------------

To build:

.. code-block:: shell

	cd spt3g_software
	mkdir build
	cd build
	cmake ..
	make


To build the documentation in the build directory type:

.. code-block:: shell

	./env-shell.sh make docs

This will construct an html version of the documentation.  This builds the documentation in the build/docs folder.  Open build/docs/index.html in your favorite web browser.  You should at least read the quick start portion of the documentation before getting started.

Release Version Tracking
------------------------

Use git tags to keep track of release versions.  Tags should be of the form "v0.1.2" for release with major version 0, minor version 1 and patch version 2.
If such a tag is defined, cmake will populate the following outputs:

 * A `cmake/Spt3gConfigVersion.cmake` file that contains the version number to be checked when including the Spt3g libraries in another cmake project
 * A `spt3g/version.py` file containing VCS parameters for access in python and stored in PipelineInfo frames
 * Add a `SPT3G_VERSION` compiler definition for accessing the version string in C++ code

When exporting the source tree to a standalone archive, run the following command in the source directory to ensure that the source version is correctly exported:

.. code-block:: shell

	cmake/config_export.sh

Then archive the source tree using  `git archive` as usual.

Version Control Hygiene
-----------------------

You can use two mechanisms to access the repository: git and SVN. The following is a brief overview of how to use these in a way that your collaborators will appreciate.

Git
===

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


SVN
===

To initially check out the repository:

.. code-block:: shell

	svn co https://user@github.com/CMB-S4/spt3g_software/trunk spt3g_software

To update your checkout:

.. code-block:: shell

	svn up

To send your changes back:

.. code-block:: shell

	svn diff files_to_commit <- Examine this
	svn ci files_to_commit

