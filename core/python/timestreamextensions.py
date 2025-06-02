import numpy
from . import G3Timestream, DoubleVector, G3VectorDouble, G3TimestreamMap, G3VectorTime, G3Time, Int64Vector, G3VectorInt, \
    G3VectorComplexDouble, ComplexDoubleVector, BoolVector, G3VectorBool, IntVector, FloatVector, ComplexFloatVector
from . import G3Units, log_fatal, log_warn, usefulfunc, G3FrameObject

__all__ = ['concatenate_timestreams']

# Use numpy bindings rather than vector_indexing_suite for G3Timestream
def G3Timestream_getitem(x, y):
    if isinstance(y, slice):
        # If the argument is a slice, and the value would be a numpy array,
        # turn it into an appropriately sliced timestream, with appropriate
        # start and stop times. This is done in C++ for speed.
        it = x._cxxslice(y)
    else:
        it = numpy.asarray(x).__getitem__(y)

    return it

@property
def timestream_dtype(ts):
    """
    Numpy data type of stored data.
    """
    return numpy.asarray(ts).dtype

G3Timestream.__getitem__ = G3Timestream_getitem
G3Timestream.__setitem__ = lambda x, y, z: numpy.asarray(x).__setitem__(y, z)
G3Timestream.dtype = timestream_dtype
G3TimestreamMap.dtype = timestream_dtype

def timestreamastype(a, dtype):
    '''
    Convert timestream to a different data type. See numpy.array.astype()
    '''
    if a.dtype == numpy.dtype(dtype):
        return a
    out = G3Timestream(numpy.asarray(a).astype(dtype))
    out.units = a.units
    out.start = a.start
    out.stop = a.stop
    return out
G3Timestream.astype = timestreamastype

def timestreammapastype(a, dtype):
    '''
    Convert timestream map to a different data type.  See numpy.array.astype()
    '''
    if a.dtype == numpy.dtype(dtype):
        return a
    return G3TimestreamMap(a.keys(), numpy.asarray(a).astype(dtype), a.start, a.stop, a.units, 0, False)
G3TimestreamMap.astype = timestreammapastype

# XXX consider replacing all of this with the numpy __array_ufunc__ machinery

# Provide all the numpy binary arithmetic operators, first non-in-place, then
# in-place, with slightly different semantics
def numpybinarywrap(a, b, op):
    is_g3 = isinstance(a, G3FrameObject) or isinstance(b, G3FrameObject)
    is_tsa = isinstance(a, G3Timestream)
    is_tsb = isinstance(b, G3Timestream)
    is_cxa = isinstance(a, (G3VectorComplexDouble, ComplexDoubleVector, ComplexFloatVector))
    is_cxb = isinstance(b, (G3VectorComplexDouble, ComplexDoubleVector, ComplexFloatVector))
    cls = a.__class__
    if is_tsa:
        if is_tsb:
            a._assert_congruence(b)
        elif is_cxb:
            return NotImplemented
    elif is_cxa:
        pass
    elif is_tsb or is_cxb:
        return NotImplemented
    out = op(numpy.asarray(a), numpy.asarray(b))
    k = out.dtype.kind
    s = out.dtype.itemsize
    if k == 'c':
        if is_g3:
            return G3VectorComplexDouble(out)
        return ComplexDoubleVector(out) if s == 16 else ComplexFloatVector(out)
    elif k == 'b':
        return G3VectorBool(out) if is_g3 else BoolVector(G3VectorBool(out))
    elif is_tsa:
        out = G3Timestream(out)
        out.units = a.units
        out.start = a.start
        out.stop = a.stop
        return out
    elif k == 'f':
        if is_g3:
            return G3VectorDouble(out)
        return DoubleVector(out) if s == 8 else FloatVector(out)
    elif k in 'iu':
        out = out.astype(int)
        if is_g3:
            return G3VectorInt(out)
        return Int64Vector(out) if s == 8 else IntVector(out)
    return NotImplemented

def numpyinplacebinarywrap(a, b, op):
    if isinstance(a, G3Timestream) and isinstance(b, G3Timestream):
        a._assert_congruence(b)
    op(numpy.asarray(a), numpy.asarray(b))
    return a

all_cls = [G3Timestream, G3VectorDouble, DoubleVector, G3VectorInt, Int64Vector,
           G3VectorComplexDouble, ComplexDoubleVector, G3VectorBool, BoolVector,
           FloatVector, IntVector, ComplexFloatVector]

for attr in ['add', 'sub', 'mul', 'div', 'truediv', 'floordiv', 'mod', 'pow',
             'and', 'or', 'xor', 'lshift', 'rshift']:
    for cls in all_cls:
        x = '__{}__'.format(attr)
        if x in numpy.ndarray.__dict__:
            setattr(cls, x,
                    lambda a, b, op=numpy.ndarray.__dict__[x]: numpybinarywrap(a, b, op))
        x = '__r{}__'.format(attr)
        if x in numpy.ndarray.__dict__:
            setattr(cls, x,
                    lambda a, b, op=numpy.ndarray.__dict__[x]: numpybinarywrap(a, b, op))
        x = '__i{}__'.format(attr)
        if x in numpy.ndarray.__dict__:
            setattr(cls, x,
                    lambda a, b, op=numpy.ndarray.__dict__[x]: numpyinplacebinarywrap(a, b, op))

# unary operators
for x in ['__neg__', '__pos__', '__invert__', '__abs__', '__bool__']:
    if x not in numpy.ndarray.__dict__:
        continue
    for cls in all_cls:
        setattr(cls, x,
                lambda a, op=numpy.ndarray.__dict__[x]: a.__class__(op(numpy.asarray(a))))

# Bind some useful nativish binary operators
for x in ['__eq__', '__ge__', '__gt__', '__le__', '__lt__', '__ne__']:
    if x not in numpy.ndarray.__dict__:
        continue
    for cls in all_cls:
        setattr(cls, x,
                lambda a, b, op=numpy.ndarray.__dict__[x]: numpybinarywrap(a, b, op))

# Bind some useful methods
for x in ["sum", "mean", "any", "all", "min", "max", "argmin", "argmax", "var", "std"]:
    if x not in numpy.ndarray.__dict__:
        continue
    for cls in all_cls:
        def ufunc_wrapper(op):
            def ufunc(a, *args, **kwargs):
                return numpy.ndarray.__dict__[op](numpy.asarray(a), *args, **kwargs)
            ufunc.__doc__ = numpy.ndarray.__dict__[op].__doc__
            return ufunc
        setattr(cls, x, ufunc_wrapper(x))

# Bind some memory-efficient methods for G3TimestreamMap
for op in ["std", "var"]:
    def ufunc_wrapper(op):
        def ufunc(a, axis=None, dtype=None, out=None, **kwargs):
            bound_args = {}
            if op in ["std", "var"]:
                bound_args["ddof"] = kwargs.pop("ddof", 0)
            if not len(kwargs) and (axis == -1 or axis == 1):
                ret = numpy.asarray(getattr(a, "_c" + op)(**bound_args))
                if out is not None:
                    out[:] = ret
                    return out
                if dtype is not None:
                    return ret.astype(dtype)
                return ret
            kwargs.update(**bound_args)
            return numpy.ndarray.__dict__[op](numpy.asarray(a), axis, dtype, out, **kwargs)
        ufunc.__doc__ = numpy.ndarray.__dict__[op].__doc__
        return ufunc

    setattr(G3TimestreamMap, op, ufunc_wrapper(op))


#add concatenation routines to g3timestream objects
def _concatenate_timestreams(cls, ts_lst, ts_rounding_error=0.6, ts_interp_threshold=0):
    """
    Concatenate G3Timestream objects together.

    Arguments
    ---------
    ts_lst : list
        list of G3Timestream objects.
    ts_rounding_error : float
        allowed error in timestream separation such that timestreams are
        contiguous, as a fraction of the sample rate. This should be
        0 by default, but is 0.5 to allow for downsampler shifting,
        and then bumpted again to 0.6 to allow for floating-point
        errors in what 0.5 is.
    ts_interp_threshold : float
        allowed timestream separation below which gaps between timestreams are
        interpolated to be made continuous

    Returns
    -------
    ts : G3Timestream instance
        The concatenation of the input list of G3Timestream objects
    """
    #check for contiguous timestreams
    for i in range(1, len(ts_lst)):
        ts_sep = (ts_lst[i].start.time - ts_lst[i-1].stop.time) * ts_lst[i].sample_rate
        if numpy.abs(ts_sep - 1) > ts_rounding_error:
            if (ts_interp_threshold > 0) and ((ts_sep - 1) < ts_interp_threshold):
                log_warn("Timestreams are not contiguous: timestreams %d and %d "
                          "separated by %f samples (%s).  Interpolating." %
                         (i, i-1, ts_sep - 1, str(ts_lst[i].start)))
                v = numpy.linspace(ts_lst[i-1][-1], ts_lst[i][0], ts_sep + 1)[1:-1]
                ts_interp = cls(v)
                ts_interp.units = ts_lst[0].units
                ts_interp.start = ts_lst[i-1].stop + 1. / ts_lst[i].sample_rate
                ts_interp.stop = ts_lst[i].start - 1. / ts_lst[i].sample_rate
                ts_lst = ts_lst[:i] + [ts_interp] + ts_lst[i:]
            else:
                log_fatal("Timestreams are not contiguous: timestreams %d and %d "
                          "separated by %f samples (%s)" % (i, i-1, ts_sep - 1, str(ts_lst[i].start)))
        if (ts_lst[i].units != ts_lst[0].units):
            log_fatal("Timestreams are not the same units")
    out_ts = cls(numpy.concatenate(ts_lst))
    out_ts.units = ts_lst[0].units
    out_ts.start = ts_lst[0].start
    out_ts.stop = ts_lst[-1].stop
    return out_ts

G3Timestream.concatenate = classmethod(_concatenate_timestreams)


def _concatenate_timestream_maps(cls, ts_map_lst, ts_rounding_error=0.6, ts_interp_threshold=0, skip_missing=False):
    """
    Concatenate G3TimestreamMap objects together.

    Arguments
    ---------
    ts_map_lst : list
        list of G3TimestreamMap objects.
    ts_rounding_error : float
        allowed error in timestream separation such that timestreams are
        contiguous, as a fraction of the sample rate. This should be
        0 by default, but is 0.5 to allow for downsampler shifting,
        and then bumpted again to 0.6 to allow for floating-point
        errors in what 0.5 is.
    ts_interp_threshold : float
        allowed timestream separation below which gaps between timestreams are
        interpolated to be made continuous
    skip_missing : bool
        If True, include only the channels that are present in all of the
        input G3TimestreamMap object.  Otherwise, raises an error if any
        map does not have the same keys.

    Returns
    -------
    tsm : G3TimestreamMap instance
        The concatenation of the input list of G3TimestreamMap objects
    """
    keys = list(ts_map_lst[0].keys())
    skeys = set(keys)
    for tsm in ts_map_lst[1:]:
        if skip_missing:
            skeys &= set(tsm.keys())
            continue
        if set(tsm.keys()) != skeys:
            log_fatal("Timestream maps do not have the same keys")
    if skip_missing:
        # preserve key order
        keys = [k for k in keys if k in skeys]
    out_tsm = cls()
    for k in keys:
        ts_lst = [tsm[k] for tsm in ts_map_lst]
        out_tsm[k] = G3Timestream.concatenate(ts_lst, ts_rounding_error, ts_interp_threshold)
    return out_tsm

G3TimestreamMap.concatenate = classmethod(_concatenate_timestream_maps)


# global concatenation function
@usefulfunc
def concatenate_timestreams(ts_lst, ts_rounding_error=0.6, ts_interp_threshold=0):
    """
    Concatenate G3Timestream or G3TimestreamMap objects together.

    Arguments
    ---------
    ts_lst : list
        list of G3Timestream or G3TimestreamMap objects.  Must all be
        the same type.
    ts_rounding_error : float
        allowed error in timestream separation such that timestreams are
        contiguous, as a fraction of the sample rate. This should be
        0 by default, but is 0.5 to allow for downsampler shifting,
        and then bumpted again to 0.6 to allow for floating-point
        errors in what 0.5 is.
    ts_interp_threshold : float
        allowed timestream separation below which gaps between timestreams are
        interpolated to be made continuous

    Returns
    -------
    ts : G3Timestream or G3TimestreamMap object
        The concatenation of the input list of objects
    """
    return ts_lst[0].concatenate(ts_lst, ts_rounding_error, ts_interp_threshold)


@property
def tsm_names(self):
    '''
    Get timestream map channel names.
    '''
    return list(self.keys())

G3TimestreamMap.names = tsm_names

@property
def tsm_data(self):
    '''
    Return a numpy array view into the underlying 2D array of the timestream map
    '''
    return numpy.asarray(self)

G3TimestreamMap.data = tsm_data
