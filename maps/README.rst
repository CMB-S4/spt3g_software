----
maps
----

The maps project defines map projections, along with G3SkyMap subclasses that provide sky maps in those map projections and tools for format/projection conversions of these data types.

In addition, this library contains three pipeline modules (MapBinner, SingleDetectorMapBinner, and SingleDetectorBoresightBinner) to make binned maps from time-ordered data, as well as a module (MapMockObserver) to mock-observe a provided sky map, generating fake time-ordered-data from it corresponding to some stored instrument pointing itinerary. A few other utility pipeline modules (described below) are provided for some map manipulation tasks.

The key data types defined here are:


HealpixSkyMap
  Implements Healpix over all or a fraction of the sky, either in nested or ring mode. The underlying sky map data are represented in one of three ways: as a dense 1-D array (full sky), as a locally-dense region surrounded by zeroes (ring mode only), or as a list of non-zero pixels and their values. The second two modes efficiently represent partial-sky maps.

FlatSkyMap
  Implements a flat-sky map (similar to a 2D numpy array) in any of the supported projections. The stored map is either a dense 2D array, or a locally dense region (a set of neighboring columns containing non-zero values, each of which contains a single contiguous block of non-zero values).

Map Attributes
==============

The following attributes are common to all G3SkyMap subclasses:

``coord_ref``
  The coordinate system on the sky to which each map pixel is referenced, stored as an instance of the ``MapCoordReference`` enum.  Currently supported coordinate systems are ``Equatorial`` (FK5 J2000), ``Galactic`` and ``Local`` (telescope azimuth and elevation).

``pol_type``
  The Stokes polarization of the map object, which is an instance of the ``MapPolType`` enum, and can have the value ``T``, ``Q``, ``U`` or None.

``pol_conv``
  The polarization convention used to encode the Q and U Stokes orientations relative to the coordinate axes.  This attribute is an instance of the ``MapPolConv`` enum, which can have the value ``IAU``, ``COSMO`` or None.  Both IAU and COSMO polarization conventions are supported in polarization-aware functions (e.g. ``FlattenPol``), but most default to using the IAU convention.  Warnings will be raised when a polarized map is used without a polarization convention set.  Changing the polarization convention between IAU and COSMO on a ``U`` map results in flipping the sign of all pixels in the map.
  
``units``
  The units system in which the map is computed, stored as an instance of the ``G3TimestreamUnits`` enum, typically ``Tcmb``.
  
``weighted``
  A boolean attribute indicating whether the data in the map have been normalized by the inverse of the appropriate Mueller matrix (``weighted=False``) or not (``weighted=True``).  See more information on map weights below.

File Format Conversions
=======================

We support writing maps into FITS files that can be read with other tools (such as DS9), using the ``fitsio.save_skymap_fits`` function.  FITS files with compatible headers can be read in using the ``fitsio.load_skymap_fits`` function.

T, Q, U and corresponding G3SkyMapWeights objects are written to a single file as a sequence of HDUs.  FlatSkyMap objects are stored in dense format (see below) to CompImageHDU objects if compression is enabled, and otherwise stored in dense format to standard ImageHDU objects.  The latter can be loaded using "old style" fits readers, such as the ``idlastro`` fits utilities.

HealpixSkyMap objects are stored in a sequence of BinTableHDU objects, a format that is compatible with the ``healpy.read_map`` function.  Dense maps (see below) are stored using implicit indexing, and sparse maps are stored using explicit indexing with an additional pixel index column.

Indexing
========

Values in maps can be set and retrieved using the standard python (or C++) ``[]`` operator. Both flat and Healpix maps support a 1-D indexing convention. For flat-sky maps, this 1-D index follows C ordering; for Healpix maps, this is the normal 1-D Healpix pixel number. Flat-sky maps also accept 2-D indices, which have ordering following normal language conventions for 2-D indices ((y, x) in Python, (x, y) in C++).

Note that sky maps *do not* support numpy-style slicing operations, except for 2-D indexing of flat-sky maps (see below), which makes a copy of the underlying map data.  To perform operations with other numpy arrays, use ``numpy.asarray``, which will convert the map to its dense representation (see below) and provides read-write access the map's internal buffer, which requires no meaningful CPU time or memory.

Sparsity
========

By default, both Healpix and flat-sky maps are initialized in sparse mode. This imposes a slight performance penalty but will result in the map storing only non-zero portions (with caveats, see details above), substantially reducing RAM usage. Some map operations, in particular casting to numpy arrays, will result in the implicit conversion of the map to dense storage, which can result in sudden increases in RAM usage. The current sparsity mode can be examined or changed with the ``sparse`` property (flat sky maps) or the ``dense``, ``ringsparse``, or ``indexedsparse`` properties (Healpix maps). Serialization to ``.g3`` files will maintain the current sparsity scheme, as do arithmetic operators where possible. Serialization to ``.fits`` files implicitly converts flat sky maps to dense mode, but preserves the sparsity of Healpix maps.  The current number of stored pixels can be obtained using the ``npix_allocated`` property, and the number of non-zero pixels can be obtained using the ``npix_nonzero`` property.  Dense maps can be efficiently compactified in memory using the ``G3SkyMap.compact`` method, or the ``CompactMaps`` pipeline module.

Beyond paying attention to implicit conversions to dense storage and the performance impact of sparse storage (which is small, at least for FlatSkyMap objects), users of this code do not need to worry about the storage mode--all interfaces are identical in all modes.

Masking
=======

Maps containing only boolean data for each pixel are stored as ``G3SkyMapMask`` objects.  Such mask objects have a ``.parent`` attribute which is a shallow clone of the map object with which they are associated (to check for shape compatibility).

Masks are returned when using comparison operators with map objects, e.g.  ``map1 > 5`` or ``map1 == map2``.  The supported comparison operators are: ``>, >=, ==, !=, <=, <``.  Masks can also be combined together using logical operators, e.g. ``mask3 = mask1 & mask2`` or ``mask1 ^= mask2``.  The supported comparison operators are: ``&, &=, |, |=, ^, ^=``.  Masks can also be checked for equality to other masks using ``==`` and ``!=`` operators.

Mask objects can be ``clone``'ed the same way as maps.  A map can be converted to a boolean mask using ``G3SkyMap.to_mask()``, which returns a mask which is ``True`` wherever the map is non-zero (optionally excluding nan or inf pixels).  A mask can be converted back to a map object using ``G3SkyMapMask.to_map()``, which returns a sparse, unit-less, unweighted, unpolarized map object of the same type as ``G3SkyMapMask.parent``, containing double ``1.0`` wherever the mask is ``True``.

Masks can also be applied to maps or masks using the appropriate ``.apply_mask`` method, with optional inversion; alternatively maps can also be directly multiplied by a compatible mask object.  A list of non-zero pixels can be returned using ``.nonzero()`` (note that this returns a single vector of pixel positions), and mask contents can be checked using ``.all()``, ``.any()`` and ``.sum()``.  Mask contents can be inverted in-place using ``.invert()``.

Mask objects cannot be accessed using ``numpy`` slicing, or converted directly to arrays, because ``numpy`` does not represent boolean values as single bits.  To be able to use ``numpy`` tools with masks, you need to first convert the mask to a dense map using ``.to_map()``.  All associated methods of the parent map are accessible as attributes of the mask object in python, e.g. ``mask.angles_to_pixels()`` works as one would expect.

Mask Memory Usage
-----------------

The current implementation of masks is to use a dense ``std::vector<bool>`` as the data storage backend, which uses 64x less memory than a dense map (``std::vector<double>``) of the same dimensions.  This implementation is sufficient for ``FlatSkyMap`` objects, since these are typically O(50\%) full populated in their sparse state; however, the memory savings for ``HealpixSkyMap`` objects is not as significant when observing sufficiently small patches of sky.  Future work would enable a similar sparse storage backend for masks.

In general, when working with high-resolution maps of any sort, it is important to think carefully about doing the sorts of operations that can balloon memory usage, e.g. taking care to preserve the sparsity of maps by avoiding numpy operations if possible, or using in-place operations to avoid unintentionally creating extra maps or masks in memory.

Statistics
==========

Most ``numpy.ufunc``-like methods are defined for map objects, namely ``.all(), .any(), .sum(), .mean(), .median(), .var(), .std(), .min(), .max(), .argmin(), .argmax()``.  All methods take an optional ``where`` argument, which can be a compatible ``G3SkyMapMask`` object, or size-compatible 1-D ``numpy`` array that can be converted into one.  In addition, these methods are called under the hood when using the numpy equivalent functions (``numpy.all()``, etc), in order to preserve the sparsity of the input map.  Methods that ignore ``NaN`` values are also defined (``.nansum()``, etc), which behave much like the standard methods, except that calling ``numpy.nansum()`` and friends on a map object does *not* preserve sparsity.

Map values can be tested using ``.isnan(), .isinf(), .isfinite()`` methods as well; these return ``G3SkyMapMask`` objects.

Map Interpolation
=================

Several interpolation and rebinning utilities are provided.  The method ``G3SkyMap.get_interp_values`` can be used for extracting map values at arbitrary sky positions using bilinear interpolation.  The method ``G3SkyMap.rebin`` can be used to downgrade the map resolution in a way that preserves the total power within each map pixel.

The functions ``healpix_to_flatsky`` and ``flatsky_to_healpix`` functions are provided to reproject maps between flat sky and curved sky systems, with options to use interpolation or rebinning to improve the accuracy of the reprojection.

The more general ``reproj_map`` function can also be used to convert between flat sky projections.

*Note:* The interpolation routine for healpix maps produces results that differ from those of the equivalent ``healpy.get_interp_val`` routine.  The interpolation routine implemented here is area-preserving in the computation of bilinear weights, whereas the ``healpy`` routine is not.

Map Weights
===========

The ``G3SkyMapWeights`` class combines the six unique components of the Mueller weight matrix into one object.  The individual matrix terms can be accessed using the attributes ``G3SkyMapWeights.TT``, etc, or as keyed elements (e.g. ``weights['TT']``).  The full matrix for an individual map pixel can be accessed using the standard ``[]`` operator.  In python, this returns a symmetric 3x3 numpy array that is a copy of the values in the underlying maps, and in C++ this returns a MuellerMatrix object, with scalar attributes ``MuellerMatrix.tt``, etc that are writable references to elements of the underlying map objects.  The ``G3SkyMapWeights.polarized`` attribute determines whether the weight structure contains polarization information.  For unpolarized weights, only the ``TT`` element is set, and the ``[]`` operator returns a scalar value in python, and a MuellerMatrix with just the TT element set in C++.

In C++ there is also a StokesVector object that is analogous to the MuellerMatrix object.  It has scalar attributes StokesVector.t etc, that are writable references to elements of map objects.  Matrix operations on the StokesVector and MuellerMatrix objects are well defined.

Weights are removed from or applied to a set of Stokes T/Q/U maps simultaneously, using the ``remove_weights`` or ``apply_weights`` functions, or their corresponding pipeline modules.

Map Frames and Pipelines
========================

Maps and associated weights are generally stored in memory and on disk in ``G3Frames`` of type ``G3FrameType.Map``, with keys ``'T', 'Q', 'U', 'Wpol'`` defined for polarized maps, and ``'T', 'Wunpol'`` defined for unpolarized maps.  Map frames can be checked for validity using the ``ValidateFrames`` pipeline module, which raises errors or warnings for missing keys or inconsistent attributes.

Map frames can be manipulated in a pipeline using some memory-efficient pipeline modules.  Weights can be applied or removed from their corresponding Stokes maps using the ``ApplyWeights`` or ``RemoveWeights`` pipeline modules.  Maps can be converted to polarized or unpolarized versions using the ``MakeMapPolarized`` and ``MakeMapUnpolarized`` modules.  They can also be compactified to their most sparse representation using the ``CompactMaps`` module.

Existing maps can be injected into a pipeline using the ``InjectMaps`` module, and map stubs can be injected using ``InjectMapStub`` or ``ReplicateMaps``.  Maps can also be extracted from a pipeline using the ``ExtractMaps`` module.

Flat Sky Map Projections
========================

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

Flat Sky Map Manipulation
=========================

Flat sky maps have additional functions defined for efficient manipulation in memory.

The ``FlattenPol`` pipeline module flattens the Q and U stokes parameters to align with the pixel coordinate grid, which is necessary for computing power spectra in the flat sky approximation.

Small patches can be extracted from and inserted into larger flat sky maps using the ``FlatSkyMap.extract_patch`` and ``FlatSkyMap.insert_patch`` methods, respectively.  Also, maps can be padded and cropped using the ``FlatSkyMap.reshape`` method, which keeps the patch centered in the output map.  All of these preserve the map pixelization and correspondence to angle on the sky.

As an equivalent and more Pythonic alternative, you can also extract portions of the map using numpy-style slicing operations (e.g. ``map[45:130,114:182]``), which will produce a map with the same contents as the numpy operation but without converting it to a dense map and with all the coordinate information set appropriately (and is equivalent to ``extract_patch()``). This also works with setting, but the coordinates have to match the sub-subcoordinates (as you would have gotten them from getting a slice or ``extract_patch()``).  Note that this slicing creates a copy of the underlying data, so in-place operations (e.g. ``map[45:130,114:182] += 5``) will work, but are not necessarily memory efficient.

Map Pointing
============

This package also provides functions and pipeline modules for creating and manipulating the quaternions necessary for mapmaking.  In general, there are two forms of quaternions that are used throughout the code: pointing quaternions and rotation quaternions.

Pointing Quaternions
--------------------

Pointing quaternions encode the two-dimensional sky coordinate angles in their vector component.  These quaternions can be created using the ``ang_to_quat`` function, and their sky coordinates extracted using the ``quat_to_ang`` function.  The various methods of the ``G3SkyMap`` classes return or accept pointing quaternions.  Note that local (horizon) coordinates have a different parity than sky coordinates (equatorial, galactic); the ``z`` vector coordinate encodes ``-sin(elevation)`` in local coordinates, but ``+sin(dec)`` in sky coordinates.

Rotation Quaternions
--------------------

Conversion between coordinate systems is done by constructing rotation quaternions.  A pointing quaternion ``q_p`` can be rotated to a new coordinate system by the rotation quaternion ``q_r`` by using quaternion multiplication: ``q_p_rot = q_r * q_p / q_r``.  For example, the module ``FillCoordTransRotations`` can be used to construct rotation quaternions for rotating detector offset coordinates into local or on-sky coordinate systems.  Rotation quaternions can be rotated into Galactic coordinates using the ``EquatorialToGalacticTransRotations`` module.

Detector Pointing
-----------------

Detector pointing timestreams are constructed by first using the ``offsets_to_quat`` function to construct the detector offset quaternion in boresight coordinates, then rotating that pointing quaternion onto the sky by applying a rotation quaternion constructed from the boresight pointing timestreams.  This is done internally for each detector in each of the mapmaking pipeline modules (``MapBinner``, ``MapMockObserver``, etc), which all require an input ``BolometerPropertiesMap`` object with offsets for each detector, and pre-computed timestreams of boresight rotation quaternions associated with each input ``Scan`` frame.

