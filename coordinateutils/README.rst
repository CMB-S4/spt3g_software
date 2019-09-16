---------------
coordinateutils
---------------

The coordinateutils project defines map projections, along with G3SkyMap subclasses that provide sky maps in those map projections and tools for format/projection conversions of these data types. 

The key data types defined here are:


HealpixSkyMap
  Implements Healpix over all or a fraction of the sky, either in nested or ring mode. The underlying sky map data are represented in one of three ways: as a dense 1-D array (full sky), as a locally-dense region surrounded by zeroes (ring mode only), or as a list of non-zero pixels and their values. The second two modes efficiently represent partial-sky maps.

FlatSkyMap
  Implements a flat-sky map (similar to a 2D numpy array) in any of the supported projections. The stored map is either a dense 2D array, or a locally dense region (a set of neighboring columns containing non-zero values, each of which contains a single contiguous block of non-zero values).

Description of Map Projections
==============================

For flat-sky maps, we support the following map projections:

ProjSansonFlamsteed
  Sanson-Flamsteed (also called the sinusoidal projection). It has equal-area pixels, defined by multiplying azimuth distances by cos(latitude). Mercator-esque in that lines of constant latitude are transformed to lines of constant y. Distances are not preserved. Also known as "proj 0".

ProjPlateCarree
  The Plate-Carree projection just plots latitude and longitude on a grid: latitude lines are at constant y and equally spaced, while longitude lines are at constant x and equally spaced. Pixels are not equal-area. Also known as "proj 1".

ProjOrthographic
  The projection of the sphere onto a plane -- the sky looks like a circle. Can only show one hemisphere. Lines drawn on the map do not correspond to latitude or longitude. Pixels are not equal-area. Also known as "proj 2".

ProjStereographic
  Another projection of the sphere onto a plane that makes it look like a circle. Differs from an orthographic projection in that it lets you see both hemispheres. Popularized in the form of the UN logo. Lines drawn on the map do not correspond to latitude or longitude. Pixels are not equal-area. Also known as "proj 4".

ProjLambertAzimuthalEqualArea
  Yet another mapping of the sphere to a circle, but this one has equal-area pixels. Largely distance-preserving, which makes it particularly useful for power-spectrum analyses. Also known as "proj 5".

ProjGnomonic
  Another projection of the sphere onto a circle. This one has the property that straight lines correspond to geodesics. Does not have equal-area pixels. Can show less than half a sphere. Also known as a "tangent projection" or "proj 6".

ProjCylindricalEqualArea
  The Lambert cylindrical equal-area projection (CEA) maps the sphere to a rectangle. Has equal-area pixels. Lines of constant x correspond to constant longitude; lines of constant y are constant latitude. Latitudes get closer together (by sin(latitude)) at the poles. Also known as "proj 7".
  

File Format Conversions
=======================

We support writing maps into FITS files that can be read with other tools (such as DS9).

Indexing
========

Values in maps can be set and retrieved using the standard python (or C++) ``[]`` operator. Both flat and Healpix maps support a 1-D indexing convention. For flat-sky maps, this 1-D index follows C ordering; for Healpix maps, this is the normal 1-D Healpix pixel number. Flat-sky maps also accept 2-D indices, which have ordering following normal language conventions for 2-D indices ((y, x) in Python, (x, y) in C++).

Note that sky maps *do not* support numpy-style slicing operations. These are ill-defined in the case of flat-sky maps (they don't preserve invariants of the projections) and of dubious utility for Healpix. If you want to take one anyway, use ``numpy.asarray`` and numpy operations. For dense maps (see below), this will access the map's internal buffer directly and so requires no meaningful CPU time or memory -- ``numpy.asarray`` also provides read/write access.

Sparsity
========

By default, both Healpix and flat-sky maps are initialized in sparse mode. This imposes a slight performance penalty but will result in the map storing only non-zero portions (with caveats, see details above), substantially reducing RAM usage. Some map operations, in particular casting to numpy arrays, will result in the implicit conversion of the map to dense storage, which can result in sudden increases in RAM usage. The current sparsity mode can be examined or changed with the ``dense`` property (flat sky maps) or the ``dense``, ``ringsparse``, or ``indexedsparse`` properties (Healpix maps). Serialization will maintain the current sparsity scheme, as do arithmetic operators where possible. The current number of stored pixels can be obtained using the ``npix_allocated`` property.

Beyond paying attention to implicit conversions to dense storage and the performance impact of sparse storage (which is small), users of this code do not need to worry about the storage mode--all interfaces are identically.

Map Interpolation
=================

Needs documentation. Apologies.

