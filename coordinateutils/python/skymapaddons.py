import numpy
from spt3g.coordinateutils import G3SkyMapWeights, G3SkyMap, WeightType

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


for x in ['__add__', '__and__', '__div__', '__divmod__', '__floordiv__', '__ge__', '__gt__', '__iadd__', '__iand__', '__idiv__', '__ifloordiv__', '__imod__', '__imul__', '__ior__', '__ipow__', '__isub__', '__itruediv__', '__le__', '__lt__', '__mul__', '__neg__', '__nonzero__', '__or__', '__pow__', '__radd__', '__rdiv__', '__rdivmod__', '__rmod__', '__rmul__', '__rpow__', '__rsub__', '__rtruediv__', '__truediv__']:
    setattr(G3SkyMap, x, lambda a, b, op=x: numpycompat(a, b, op))

# Bind all numpy unary operators to G3SkyMap
for x in ['__neg__', '__pos__', '__invert__']:
    setattr(G3SkyMap, x, lambda a, op=x: numpyinplace(a, numpy.ndarray.__dict__[op](numpy.asarray(a))))

# And some special binary operators
setattr(G3SkyMap, '__getslice__', lambda a, *args: numpy.ndarray.__getslice__(numpy.asarray(a), *args))

# Make weight maps so that you can index them and get the full 3x3 weight matrix

def skymapweights_getitem(self, x):
    mat = numpy.zeros((3,3))
    mat[0,0] = self.TT[x]
    mat[0,1] = mat[1,0] = self.TQ[x]
    mat[0,2] = mat[2,0] = self.TU[x]
    mat[1,1] = self.QQ[x]
    mat[1,2] = mat[2,1] = self.QU[x]
    mat[2,2] = self.QU[x]

    return mat

G3SkyMapWeights.__getitem__ = skymapweights_getitem
del skymapweights_getitem

def skymapweights_setitem(self, x, mat):
    # Check for symmetry, shape?
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
        setattr(self.TQ, x, val)
        setattr(self.TU, x, val)
        setattr(self.QQ, x, val)
        setattr(self.QU, x, val)
        setattr(self.UU, x, val)
G3SkyMapWeights.__setattr__ = skymapweights_setattr
del skymapweights_setattr

# Bind addition to weights
def add_skymapweights(w1,w2):
    w3 = G3SkyMapWeights(w1.TT,w1.weight_type)
    w3.TT = w1.TT + w2.TT

    if w1.weight_type == WeightType.Wpol:
        w3.TQ = w1.TQ + w2.TQ
        w3.TU = w1.TU + w2.TU
        w3.QQ = w1.QQ + w2.QQ
        w3.QU = w1.QU + w2.QU
        w3.UU = w1.UU + w2.UU

    return w3

G3SkyMapWeights.__add__ = add_skymapweights
del add_skymapweights

