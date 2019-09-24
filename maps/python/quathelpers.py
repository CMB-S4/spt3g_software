from spt3g.core import G3Units, quat, G3VectorQuat
import numpy

def quat_to_ang(q):
    '''
    Convert a pointing quaternion (or vector of them) to a set of angles
    (or vector of them) specified as a (longitude, latitude) pair.
    '''
    single = False
    if isinstance(q, quat):
        q = numpy.asarray(G3VectorQuat([q]))
        single = True
    elif isinstance(q, list):
        q = numpy.asarray(G3VectorQuat(q))
    else:
        q = numpy.asarray(q)

    # Copied from C code
    d = (q[:,1]**2 + q[:,2]**2 + q[:,3]**2)
    mask = numpy.abs(d - 1.0) > 1e-6
    q[mask] /= d[mask][:,None]**.5

    delta = numpy.arcsin(q[:,3])*G3Units.rad
    alpha = numpy.arctan2(q[:,2], q[:,1])*G3Units.rad

    if single:
        return (alpha[0], delta[0])
    else:
        return (alpha, delta)

def ang_to_quat(alpha, delta):
    '''
    Convert a set of angles (or vector of them) specified as a
     (longitude, latitude) pair to a pointing quaternion (or vector of them).
    '''

    alpha = numpy.asarray(alpha)
    delta = numpy.asarray(delta)
    # Copied from C code
    c_delta = numpy.cos(delta/G3Units.rad)
    q = numpy.column_stack((0*c_delta, # 0s with the right shape
        c_delta * numpy.cos(alpha/G3Units.rad),
        c_delta * numpy.sin(alpha/G3Units.rad),
        numpy.sin(delta/G3Units.rad)))

    if len(q) == 1:
        return quat(q[0,0], q[0,1], q[0,2], q[0,3])
    else:
        return G3VectorQuat(q)

