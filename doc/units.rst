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

Time:

* nanosecond, ns
* microsecond, us
* millisecond, ms
* second, s
* minute, min
* hour, h
* day

Frequency:

* Hz, hz
* MHz
* GHz

Angles:

* rad
* deg
* arcmin
* arcsec
* rahour, rahr

Length:

* nanometer, nm
* micron
* millimeter, mm
* centimeter, cm
* inch, in
* foot, ft
* meter, m
* kilometer, m
* au
* parsec, pc

Power:

* attowatt, aW
* picowatt, pW
* nanowatt, nW
* microwatt, uW
* milliwatt, mW
* watt, W
* horsepower, hp

Voltage:

* volt, V
* millivolt, mV
* microvolt, uV

Current:

* amp, ampere, A
* milliamp, mA
* microamp, uA
* nanoamp, nA

Temperature:

* picokelvin, pK
* nanokelvin, nK
* microkelvin, uK
* millikelvin, mK
* kelvin, K
* rankine, R

