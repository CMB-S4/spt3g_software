from ..core import G3Units, Quat, G3VectorQuat, G3TimestreamQuat, usefulfunc, indexmod
import numpy


@usefulfunc
def quat_to_ang(q):
    """
    Convert a pointing quaternion (or vector of them) to a set of angles (or
    vector of them) specified as a (longitude, latitude) pair.
    """
    single = False
    if isinstance(q, Quat):
        q = numpy.array(G3VectorQuat([q]))
        single = True
    elif isinstance(q, list):
        q = numpy.array(G3VectorQuat(q))
    else:
        q = numpy.array(q)

    # Copied from C code
    d = q[:, 1] ** 2 + q[:, 2] ** 2 + q[:, 3] ** 2
    mask = numpy.abs(d - 1.0) > 1e-6
    q[mask] /= d[mask][:, None] ** 0.5

    delta = numpy.arcsin(q[:, 3]) * G3Units.rad
    alpha = numpy.arctan2(q[:, 2], q[:, 1]) * G3Units.rad

    if single:
        return (alpha[0], delta[0])
    else:
        return (alpha, delta)


@usefulfunc
def ang_to_quat(alpha, delta, start=None, stop=None):
    """
    Convert a set of angles (or vector of them) specified as a (longitude,
    latitude) pair to a pointing quaternion (or vector of them). If start
    and stop are defined, the return value for vectors will be a
    G3TimestreamQuat with start and stop times set to the provided values.
    """

    alpha = numpy.asarray(alpha) / G3Units.rad
    delta = numpy.asarray(delta) / G3Units.rad
    # Copied from C code
    c_delta = numpy.cos(delta)
    q = numpy.column_stack(
        (
            0 * c_delta,  # 0s with the right shape
            c_delta * numpy.cos(alpha),
            c_delta * numpy.sin(alpha),
            numpy.sin(delta),
        )
    )

    if len(q) == 1:
        return Quat(q[0, :])
    else:
        if start is not None:
            out = G3TimestreamQuat(q)
            out.start = start
            out.stop = stop
            return out
        else:
            return G3VectorQuat(q)

@indexmod
def AddTimingToPointingQuats(fr, key, timing_ref='RawBoresightAz'):
    """
    Transforms a quaternion pointing timestream <key> from a G3VectorQuat
    to a G3TimestreamQuat by adding timing information extracted from some
    other timestream-like object (scalar pointing timestreams, detector
    timestreams etc.) specified by <timing_ref>. Because this operation is
    backwards-compatible and involves no data loss, the transformation is
    done in-place -- the previous data are deleted and replaced with the new
    one.
    """
    
    if key not in fr:
        return

    if isinstance(fr[key], G3TimestreamQuat):
        return

    x = fr.pop(key)
    x = G3TimestreamQuat(x)
    x.start = fr[timing_ref].start
    x.stop = fr[timing_ref].stop
    fr[key] = x

