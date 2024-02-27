import numpy
import warnings
from spt3g.maps import G3SkyMapWeights, G3SkyMap, FlatSkyMap, G3SkyMapMask

# This file adds extra functionality to the python interface to G3SkyMap and
# G3SkyMapWeights. This is done in ways that exploit a large fraction of
# the evil things you can do in Python (Nick referred to it as "devil magic"),
# but please believe that the alternative is worse. I have tried to explain
# the rationale for each block.

# Bind all numpy binary operators to G3SkyMap by stealing the numpy functions
# and wrapping them in short lambdas that cast to and from sky maps

def numpyinplace(skymap, array):
    # Some magic to give skymaps the same operators as numpy arrays. This
    # requires multiple lines, so cannot be in the lambda below.  Because
    # this calls __copy__(), which is Clone(), it does the right thing for
    # subclasses of skymaps.
    rv = skymap.__copy__()
    numpy.asarray(rv)[:] = array
    return rv

def numpycompat(a, b, op):
    if isinstance(b, G3SkyMap) and not a.compatible(b):
        raise TypeError("Map of type <{}> is incompatible with map of type <{}>".format(a, b))
    return numpyinplace(a, numpy.ndarray.__dict__[op](numpy.asarray(a), numpy.asarray(b)))


for x in ['__and__', '__divmod__', '__floordiv__', '__iand__', '__ifloordiv__', '__imod__', '__ior__', '__mod__', '__or__', '__rdivmod__', '__rmod__', '__rpow__']:
    setattr(G3SkyMap, x, lambda a, b, op=x: numpycompat(a, b, op))

for op in ["all", "any", "sum", "mean", "var", "std", "min", "max", "argmin", "argmax"]:
    def ufunc_wrapper(op):
        def ufunc(a, *args, **kwargs):
            bound_args = {}
            bound = getattr(a, "_c" + op)
            if op in ["std", "var", "nanstd", "nanvar"]:
                if "ddof" in kwargs:
                    bound_args["ddof"] = kwargs.pop("ddof")
            if "where" in kwargs:
                where = kwargs.pop("where")
                if isinstance(where, numpy.ndarray):
                    where = G3SkyMapMask(a, where)
                bound_args["where"] = where
            for k in ["axis", "out", "dtype"]:
                if k in kwargs and kwargs[k] is None:
                    kwargs.pop(k)
            if len(args) or len(kwargs):
                raise TypeError("ufunc {} does not support complex arguments".format(op))
            return bound(**bound_args)
        ufunc.__doc__ = getattr(numpy.ndarray, op).__doc__
        return ufunc
    setattr(G3SkyMap, op, ufunc_wrapper(op))
    if op in ["all", "any", "sum"]:
        setattr(G3SkyMapMask, op, ufunc_wrapper(op))

# Make weight maps so that you can index them and get the full 3x3 weight matrix

def skymapweights_keys(self):
    """
    Return the list of string names of valid weight attributes.  For unpolarized
    weights, this list includes only TT.  Otherwise, the list includes all six
    unique weight attributes in row major order: TT, TQ, TU, QQ, QU, UU.
    """
    if self.polarized:
        return ['TT', 'TQ', 'TU', 'QQ', 'QU', 'UU']
    return ['TT']

G3SkyMapWeights.keys = skymapweights_keys
del skymapweights_keys

def skymapweights_getitem(self, x):
    if isinstance(x, str) and x in self.keys():
        return getattr(self, x)

    if (isinstance(x, tuple) and any(isinstance(xx, slice) for xx in x)) \
       or isinstance(x, slice):
        out = self.__class__()
        out.TT = self.TT[x]
        if self.polarized:
            out.TQ = self.TQ[x]
            out.TU = self.TU[x]
            out.QQ = self.QQ[x]
            out.QU = self.QU[x]
            out.UU = self.UU[x]
        return out

    if not self.polarized:
        return self.TT[x]

    mat = numpy.zeros((3,3))
    mat[0,0] = self.TT[x]
    mat[0,1] = mat[1,0] = self.TQ[x]
    mat[0,2] = mat[2,0] = self.TU[x]
    mat[1,1] = self.QQ[x]
    mat[1,2] = mat[2,1] = self.QU[x]
    mat[2,2] = self.UU[x]

    return mat

G3SkyMapWeights.__getitem__ = skymapweights_getitem
del skymapweights_getitem

def skymapweights_setitem(self, x, mat):
    if isinstance(x, str) and x in self.keys():
        setattr(self, x, mat)
        return

    if (isinstance(x, tuple) and any(isinstance(xx, slice) for xx in x)) \
       or isinstance(x, slice):
        raise NotImplementedError

    if not self.polarized:
        assert(numpy.isscalar(mat))
        self.TT[x] = mat
        return

    mat = numpy.asarray(mat, dtype=float)
    assert(mat.shape == (3, 3))
    assert(numpy.allclose(mat, mat.T))
    self.TT[x] = mat[0,0]
    self.TQ[x] = mat[0,1]
    self.TU[x] = mat[0,2]
    self.QQ[x] = mat[1,1]
    self.QU[x] = mat[1,2]
    self.UU[x] = mat[2,2]

G3SkyMapWeights.__setitem__ = skymapweights_setitem
del skymapweights_setitem

# patch handling for flat sky maps
def skymapweights_extract(self, x0, y0, width, height, fill=0):
    if not isinstance(self.TT, FlatSkyMap):
        raise NotImplementedError

    w = self.__class__()
    for k in self.keys():
        setattr(w, k, self[k].extract_patch(x0, y0, width, height, fill=fill))

    return w

skymapweights_extract.__doc__ = FlatSkyMap.extract_patch.__doc__
G3SkyMapWeights.extract_patch = skymapweights_extract
del skymapweights_extract

def skymapweights_insert(self, weights_patch, ignore_zeros=False):
    if not isinstance(self.TT, FlatSkyMap):
        raise NotImplementedError

    assert weights_patch.polarized == self.polarized, \
        "Input weights patch must be {}polarized".format("" if self.polarized else "un")

    for k in self.keys():
        self[k].insert_patch(weights_patch[k], ignore_zeros=ignore_zeros)

skymapweights_insert.__doc__ = FlatSkyMap.insert_patch.__doc__
G3SkyMapWeights.insert_patch = skymapweights_insert
del skymapweights_insert

def skymapweights_reshape(self, width, height, fill=0):
    if not isinstance(self.TT, FlatSkyMap):
        raise NotImplementedError

    w = self.__class__()
    for k in self.keys():
        setattr(w, k, self[k].reshape(width, height, fill=fill))

    return w

skymapweights_reshape.__doc__ = FlatSkyMap.reshape.__doc__
G3SkyMapWeights.reshape = skymapweights_reshape
del skymapweights_reshape

# Pass through attributes to submaps. This is not ideal because the properties
# are hidden and not visible by tab completion. A better solution would use
# __getattribute__
def skymapweights_getattr(self, x):
    if x == "flat_pol" and isinstance(self.TQ, FlatSkyMap):
        return getattr(self.TQ, x)
    if x == "pol_conv" and isinstance(self.TU, FlatSkyMap):
        return getattr(self.TU, x)
    if hasattr(self.TT, x):
        return getattr(self.TT, x)
    raise AttributeError("'{}' object has no attribute '{}'".format(type(self), x))
G3SkyMapWeights.__getattr__ = skymapweights_getattr
del skymapweights_getattr
oldsetattr = G3SkyMapWeights.__setattr__
def skymapweights_setattr(self, x, val):
    try:
        oldsetattr(self, x, val)
    except AttributeError:
        if self.TT is not None:
            setattr(self.TT, x, val)
        if self.TQ is not None:
            setattr(self.TQ, x, val)
        if self.TU is not None:
            setattr(self.TU, x, val)
        if self.QQ is not None:
            setattr(self.QQ, x, val)
        if self.QU is not None:
            setattr(self.QU, x, val)
        if self.UU is not None:
            setattr(self.UU, x, val)
G3SkyMapWeights.__setattr__ = skymapweights_setattr
del skymapweights_setattr

def skymapmask_getattr(self, x):
    if hasattr(self.parent, x):
        return getattr(self.parent, x)
    raise AttributeError("'{}' object has no attribute '{}'".format(type(self), x))
G3SkyMapMask.__getattr__ = skymapmask_getattr
del skymapmask_getattr

def skymap_clone(self, copy_data=True):
    with warnings.catch_warnings():
        warnings.simplefilter("default")
        warnings.warn(
            "{0}.Clone is deprecated, use {0}.clone instead".format(
                self.__class__.__name__
            ),
            DeprecationWarning,
        )
    return self.clone(copy_data)
G3SkyMap.Clone = skymap_clone
G3SkyMapWeights.Clone = skymap_clone
del skymap_clone

def skymap_compat(self, other):
    with warnings.catch_warnings():
        warnings.simplefilter("default")
        warnings.warn(
            "{0}.IsCompatible is deprecated, use {0}.compatible instead".format(
                self.__class__.__name__
            ),
            DeprecationWarning,
        )
    return self.compatible(other)
G3SkyMap.IsCompatible = skymap_compat
del skymap_compat
