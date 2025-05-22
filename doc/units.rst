Units
-----

This software includes a units system that is meant to end wondering whether a given function takes radians or degrees as an argument, or whether a stored time is in milliseconds or seconds. The support code is accessible to both C++ and Python as part of the ``G3Units`` namespace (``core.G3Units.X`` in Python and ``G3Units::X`` in C++).

Defining a quantity with units
==============================

Defining a quantity with units is meant to work the way you would write it down. The period of time 5 seconds would be represented in Python as:

.. code-block:: python

    5*core.G3Units.s

Under the hood, these are multiplicative constants that have the appropriate scalings relative to one another (e.g. ``core.G3Units.s`` is a thousand times larger than ``core.G3Units.ms``).

Converting to different units
=============================

If you want a numeric quantity in specific units, you divide by the unit you want (think of doing unit analysis in freshman physics courses).

For example, to take a passed angle and return its sin, irrespective of the user's preferred units:

.. code-block:: python

    def sinofstuff(angle):
        return math.sin(angle/core.G3Units.rad)

    print('sin(5 degrees): %f' % sinoffstuff(5*core.G3Units.deg))

Available units
===============

We currently define units for times, angles, lengths, power, voltage, current, and temperature. Some of these have aliases with abbreviated names (e.g. both ms and millisecond). Where this occurs, all names for a unit are separated with commas. For a value ``X`` in the following lists, it is exposed to Python as ``core.G3Units.X`` and to C++ as ``G3Units::X``.

Time
~~~~

* nanosecond(s), ns
* microsecond(s), us
* millisecond(s), ms
* second(s), sec, s
* minute(s), min
* hour(s), h
* day(s)

Frequency
~~~~~~~~~

* Hz, hz
* kHz
* MHz
* GHz

Angles
~~~~~~

* radian(s), rad
* degree(s), deg
* arcmin
* arcsec
* rahour, rahr
* raminute
* rasecond

Length
~~~~~~

* nanometer, nm
* micron
* millimeter, mm
* centimeter, cm
* inch, in
* foot, ft
* meter(s), m
* kilometer, m
* AU, au
* parsec, pc

Power
~~~~~

* attowatt, aW
* picowatt, pW
* nanowatt, nW
* microwatt, uW
* milliwatt, mW
* watt, W
* horsepower, hp

Flux Density
~~~~~~~~~~~~

* jansky, Jy
* millijansky, mJy
* megajansky, MJy

Solid Angle
~~~~~~~~~~~

* steradian(s), sr
* squaredegree(s), sqdeg, deg2
* squarearcmin, sqarcmin, arcmin2

Voltage
~~~~~~~

* volt, V
* millivolt, mV
* microvolt, uV

Current
~~~~~~~

* amp, ampere, A
* milliamp, mA
* microamp, uA
* nanoamp, nA

Resistance
~~~~~~~~~~

* ohm
* milliohm, mohm

Temperature
~~~~~~~~~~~

* picokelvin, pK
* nanokelvin, nK
* microkelvin, uK
* millikelvin, mK
* kelvin, K
* rankine, R

Pressure
~~~~~~~~

* bar, b
* millibar, mb
* Pa
* kPa

Data Size
~~~~~~~~~

* bit
* byte, B
* kilobyte, KB
* megabyte, MB
* gigabyte, GB

Mass
~~~~

* gram, g
* kilogram, kg
* milligram, mg


Common scientific constants
===========================

Some common constants are defined and made available in the ``core.G3Constants`` namespace.  These are multiplied by the appropriate combinations of the above G3Units.

Available constants are:

* ``c``: speed of light
* ``h``, ``hbar``: Planck constant and reduced Planck constant
* ``k``, ``kb``: Boltzmann constant
* ``G``: Newtonian gravitational constant
* ``g0``: Standard gravitational acceleration
* ``e``: Elementary charge
