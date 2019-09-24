----
maps
----

The maps project defines map projections, along with G3SkyMap subclasses that provide sky maps in those map projections and tools for format/projection conversions of these data types. 

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
  The Plate-Carree projection just plots latitude and longitude on a grid: latitude lines are at constant y and equally spaced, while longitude lines are at constant x and equally spaced. Pixels are not equal-area. Also known as "proj 1".  A variant of this projection, called ProjBICEP (or "proj 9"), adjusts the resolution along x to scale with the cosine of the latitude of the center of the map.

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

We support writing maps into FITS files that can be read with other tools (such as DS9), using the ``fitsio.save_skymap_fits`` function.  FITS files with compatible headers can be read in using the ``fitsio.load_skymap_fits`` function.

T, Q, U and corresponding G3SkyMapWeights objects are written to a single file as a sequence of HDUs.  FlatSkyMap objects are stored in dense format (see below) to CompImageHDU objects if compression is enabled, and otherwise stored in dense format to standard ImageHDU objects.  The latter can be loaded using "old style" fits readers, such as the ``idlastro`` fits utilities.

HealpixSkyMap objects are stored in a sequence of BinTableHDU objects, a format that is compatible with the ``healpy.read_map`` function.  Dense maps (see below) are stored using implicit indexing, and sparse maps are stored using explicit indexing with an additional pixel index column.

Indexing
========

Values in maps can be set and retrieved using the standard python (or C++) ``[]`` operator. Both flat and Healpix maps support a 1-D indexing convention. For flat-sky maps, this 1-D index follows C ordering; for Healpix maps, this is the normal 1-D Healpix pixel number. Flat-sky maps also accept 2-D indices, which have ordering following normal language conventions for 2-D indices ((y, x) in Python, (x, y) in C++).

Note that sky maps *do not* support numpy-style slicing operations. These are ill-defined in the case of flat-sky maps (they don't preserve invariants of the projections) and of dubious utility for Healpix. If you want to take one anyway, use ``numpy.asarray`` and numpy operations. For dense maps (see below), this will access the map's internal buffer directly and so requires no meaningful CPU time or memory -- ``numpy.asarray`` also provides read/write access.

Sparsity
========

By default, both Healpix and flat-sky maps are initialized in sparse mode. This imposes a slight performance penalty but will result in the map storing only non-zero portions (with caveats, see details above), substantially reducing RAM usage. Some map operations, in particular casting to numpy arrays, will result in the implicit conversion of the map to dense storage, which can result in sudden increases in RAM usage. The current sparsity mode can be examined or changed with the ``sparse`` property (flat sky maps) or the ``dense``, ``ringsparse``, or ``indexedsparse`` properties (Healpix maps). Serialization to ``.g3`` files will maintain the current sparsity scheme, as do arithmetic operators where possible. Serialization to ``.fits`` files implicitly converts flat sky maps to dense mode, but preserves the sparsity of Healpix maps.  The current number of stored pixels can be obtained using the ``npix_allocated`` property.

Beyond paying attention to implicit conversions to dense storage and the performance impact of sparse storage (which is small), users of this code do not need to worry about the storage mode--all interfaces are identical in all modes.

Map Interpolation
=================

Several interpolation and rebinning utilities are provided.  The method ``G3SkyMap.get_interp_values`` can be used for extracting map values at arbitrary sky positions.  The method ``G3SkyMap.rebin`` can be used to downgrade the map resolution in a way that preserves the total power within each map pixel.

The functions ``maputils.healpix_to_flatsky`` and ``maputils.flatsky_to_healpix`` functions are provided to reproject maps between flat sky and curved sky systems, with options to use interpolation or rebinning to improve the accuracy of the reprojection.

The more general ``maputils.reproj`` function can also be used to convert between flat sky projections.

Stokes Vectors and Mueller Matrices
===================================

We provide two additional map classes that combine T/Q/U Stokes maps and Mueller weight matrices into a single structure.

The G3SkyMapWeights class combines the six unique components of the Mueller weight matrix into one object.  The individual matrix terms can be accessed using the attributes G3SkyMapWeights.TT, etc.  The full matrix for an individual pixel can be accessed using the standard ``[]`` operator.  In python, this returns a symmetric 3x3 numpy array that is a copy of the values in the underlying maps, and in C++ this returns a MuellerMatrix object, with scalar attributes MuellerMatrix.tt, etc that are writable references to elements of the underlying map objects.

Analogously, the G3SkyMapWithWeights class combines the T, Q and U Stokes maps and the corresponding weight matrix into a single structure.  The individual Stokes components can be accessed using the attributes G3SkyMapWithWeights.T etc, and the 3-element StokesVector can be accessed with the standard ``[]`` operator.  In python, this returns a numpy array of length 3 that is a copy of the values in the underlying maps, and in C++ this returns a StokesVector object with scalar attributes StokesVector.t etc, that are writable references to elements of the underlying map objects.  Weights can be applied to or removed from the Stokes maps using the ``G3SkyMapWithWeights.apply_weights`` and ``G3SkyMapWithWeights.remove_weights`` methods.  The attribute ``G3SkyMapWithWeights.weighted`` is True if the Stokes components are weighted, and the attribute ``G3SkyMapWithWeights.polarized`` is True if the structure contains Q and U Stokes components.

Several convenience functions and operators are provided for the StokesVector and MuellerMatrix objects in C++, including typical matrix operations that are performed with these two structures.  A limited set of operators are also provided in both C++ and python for the G3SkyMapWeights and G3SkyMapWithWeights objects (addition of two such objects [e.g. for coadding], multiplication/division by a scalar, and multiplication by a G3SkyMap object [e.g. for masking]).
