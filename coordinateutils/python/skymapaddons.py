import numpy
from spt3g.coordinateutils import G3SkyMapWeights, G3SkyMapWithWeights, G3SkyMap, WeightType

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
    if isinstance(b, G3SkyMap) and not a.IsCompatible(b):
        raise TypeError("Map of type <{}> is incompatible with map of type <{}>".format(a, b))
    return numpyinplace(a, numpy.ndarray.__dict__[op](numpy.asarray(a), numpy.asarray(b)))


for x in ['__and__', '__divmod__', '__floordiv__', '__ge__', '__gt__', '__iand__', '__ifloordiv__', '__imod__', '__ior__', '__ipow__', '__le__', '__lt__', '__nonzero__', '__or__', '__pow__', '__rdivmod__', '__rmod__', '__rpow__']:
    setattr(G3SkyMap, x, lambda a, b, op=x: numpycompat(a, b, op))

# Bind all numpy unary operators to G3SkyMap
for x in ['__pos__', '__invert__']:
    setattr(G3SkyMap, x, lambda a, op=x: numpyinplace(a, numpy.ndarray.__dict__[op](numpy.asarray(a))))

# And some special binary operators
setattr(G3SkyMap, '__getslice__', lambda a, *args: numpy.ndarray.__getslice__(numpy.asarray(a), *args))

# Make weight maps so that you can index them and get the full 3x3 weight matrix

def skymapweights_getitem(self, x):
    if self.weight_type == WeightType.Wunpol:
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
    if self.weight_type == WeightType.Wunpol:
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

# Pass through attributes to submaps. This is not ideal because the properties
# are hidden and not visible by tab completion. A better solution would use
# __getattribute__
G3SkyMapWeights.__getattr__ = lambda self, x: getattr(self.TT, x)
oldsetattr = G3SkyMapWeights.__setattr__
def skymapweights_setattr(self, x, val):
    try:
        oldsetattr(self, x, val)
    except AttributeError:
        setattr(self.TT, x, val)
        if self.weight_type == WeightType.Wpol:
            setattr(self.TQ, x, val)
            setattr(self.TU, x, val)
            setattr(self.QQ, x, val)
            setattr(self.QU, x, val)
            setattr(self.UU, x, val)
G3SkyMapWeights.__setattr__ = skymapweights_setattr
del skymapweights_setattr

# Make maps with weights so that you can index them and get the full 1x3 Stokes vector

def skymapwithweights_getitem(self, x):
    if not self.polarized:
        return self.T[x]

    vec = numpy.zeros(3)
    vec[0] = self.T[x]
    vec[1] = self.Q[x]
    vec[2] = self.U[x]
    return vec

G3SkyMapWithWeights.__getitem__ = skymapwithweights_getitem
del skymapwithweights_getitem

def skymapwithweights_setitem(self, x, vec):
    if not self.polarized:
        assert(numpy.isscalar(vec))
        self.T[x] = vec
        return

    vec = numpy.asarray(vec, dtype=float)
    assert(len(vec) == 3)
    self.T[x] = vec[0]
    self.Q[x] = vec[1]
    self.U[x] = vec[2]

G3SkyMapWithWeights.__setitem__ = skymapwithweights_setitem
del skymapwithweights_setitem

# Pass through attributes to submaps. This is not ideal because the properties
# are hidden and not visible by tab completion. A better solution would use
# __getattribute__
G3SkyMapWithWeights.__getattr__ = lambda self, x: getattr(self.T, x)
oldsetattr = G3SkyMapWithWeights.__setattr__
def skymapwithweights_setattr(self, x, val):
    try:
        oldsetattr(self, x, val)
    except AttributeError:
        setattr(self.T, x, val)
        if self.polarized:
            setattr(self.Q, x, val)
            setattr(self.U, x, val)
        if self.weighted:
            setattr(self.weights, x, val)
G3SkyMapWithWeights.__setattr__ = skymapwithweights_setattr
del skymapwithweights_setattr

