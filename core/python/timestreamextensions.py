import numpy
from spt3g.core import G3Timestream, DoubleVector, G3VectorDouble, G3TimestreamMap, G3Time, usefulfunc
from spt3g.core import G3Units, log_fatal, log_warn

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
def timestream_shape(ts):
    return numpy.asarray(ts).shape

@property
def timestream_ndim(ts):
    return len(numpy.asarray(ts).shape)

G3Timestream.__getitem__ = G3Timestream_getitem
G3Timestream.__setitem__ = lambda x, y, z: numpy.asarray(x).__setitem__(y, z)
G3Timestream.__len__ = lambda x: numpy.asarray(x).__len__()
G3Timestream.shape = timestream_shape
G3Timestream.ndim = timestream_ndim

# Provide all the numpy binary arithmetic operators, first non-in-place, then
# in-place, with slightly different semantics
def numpybinarywrap(a, b, op):
    if isinstance(b, G3Timestream):
        a._assert_congruence(b)
    out = G3Timestream(op(numpy.asarray(a), numpy.asarray(b)))
    out.units = a.units
    out.start = a.start
    out.stop = a.stop
    return out

for x in ['__add__', '__and__', '__div__', '__divmod__', '__floordiv__', 
          '__mul__', '__neg__', '__or__', '__pow__', '__sub__', '__radd__',
          '__rdiv__', '__rdivmod__', '__rmod__', '__rmul__', '__rpow__',
          '__rsub__', '__rtruediv__', '__truediv__']:
    if x in numpy.ndarray.__dict__:
        setattr(G3Timestream, x, 
                lambda a, b, op=numpy.ndarray.__dict__[x]: numpybinarywrap(a, b, op))

def numpyinplacebinarywrap(a, b, op):
    if isinstance(b, G3Timestream):
        a._assert_congruence(b)
    op(numpy.asarray(a), numpy.asarray(b))
    return a

for x in ['__iadd__', '__iand__', '__idiv__', '__ifloordiv__', '__imod__', 
          '__imul__', '__ior__', '__ipow__', '__isub__', '__itruediv__']:
    if x in numpy.ndarray.__dict__:
        setattr(G3Timestream, x, 
                lambda a, b, op=numpy.ndarray.__dict__[x]: numpyinplacebinarywrap(a, b, op))

# Bind some useful nativish binary operators
for x in ['__eq__', '__ge__', '__gt__', '__le__', '__lt__', '__neq__']:
    if x in numpy.ndarray.__dict__:
        setattr(G3Timestream, x, 
                lambda a, b, op=x: numpy.ndarray.__dict__[op](numpy.asarray(a), numpy.asarray(b)))
G3Timestream.mean = lambda a, *args, **kwargs: numpy.ndarray.mean(numpy.asarray(a), *args, **kwargs)
G3Timestream.min = lambda a, *args, **kwargs: numpy.ndarray.min(numpy.asarray(a), *args, **kwargs)
G3Timestream.max = lambda a, *args, **kwargs: numpy.ndarray.max(numpy.asarray(a), *args, **kwargs)

# Bind a few useful function to double vectors (more or less timestreams) so
# that unit conversions are easy.
DoubleVector.__sub__ = lambda self, val: DoubleVector(numpy.asarray(self)-val)
DoubleVector.__add__ = lambda self, val: DoubleVector(numpy.asarray(self)+val)
DoubleVector.__mul__ = lambda self, val: DoubleVector(numpy.asarray(self)*val)
DoubleVector.__div__ = lambda self, val: DoubleVector(numpy.asarray(self)/val)
DoubleVector.__truediv__ = DoubleVector.__div__

# Prevent G3VectorDoubles from turning into DoubleVector on mult. C++ has no
# idea what to do with them.
G3VectorDouble.__add__ = lambda self, val: G3VectorDouble(numpy.asarray(self)+val)
G3VectorDouble.__sub__ = lambda self, val: G3VectorDouble(numpy.asarray(self)-val)
G3VectorDouble.__mul__ = lambda self, val: G3VectorDouble(numpy.asarray(self)*val)
G3VectorDouble.__div__ = lambda self, val: G3VectorDouble(numpy.asarray(self)/val)
G3VectorDouble.__truediv__ = G3VectorDouble.__div__

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


def _concatenate_timestream_maps(cls, ts_map_lst, ts_rounding_error=0.6, ts_interp_threshold=0):
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

    Returns
    -------
    tsm : G3TimestreamMap instance
        The concatenation of the input list of G3TimestreamMap objects
    """

    for tsm in ts_map_lst[1:]:
        if set(tsm.keys()) != set(ts_map_lst[0].keys()):
            log_fatal("Timestream maps do not have the same keys")
    out_tsm = cls()
    for k in ts_map_lst[0].keys():
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


def timestream_elapsed(self):
    '''
    Compute elapsed time array for samples.
    '''
    if self.n_samples == 1:
        return numpy.array([0]).astype(int)
    return (numpy.arange(self.n_samples) / self.sample_rate).astype(int)

G3Timestream.elapsed = timestream_elapsed
G3TimestreamMap.elapsed = timestream_elapsed


def timestream_t(self):
    '''
    Compute time vector for samples.
    '''
    times = self.elapsed() + self.start.time
    return [G3Time(t) for t in times]

G3Timestream.times = timestream_t
G3TimestreamMap.times = timestream_t

