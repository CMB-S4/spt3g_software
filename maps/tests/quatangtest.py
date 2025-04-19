#!/usr/bin/env python
from spt3g import core, maps
import numpy as np

a = core.Quat(2,3,4,5)
a2 = core.Quat(2,1,8,5)
b = core.G3VectorQuat([a, a**2, 2*a, a2])

assert(np.allclose(maps.quat_to_ang(a), maps.c_quat_to_ang_(a)))
assert(np.allclose(maps.quat_to_ang(b), np.asarray([maps.c_quat_to_ang_(x) for x in b]).transpose()))

angle = (.3, 0.4)
angles = ((.3, 0.8), (0.4, 0.2))

q = maps.ang_to_quat(angle[0], angle[1])
p = maps.c_ang_to_quat_(angle[0], angle[1])
assert(np.allclose(core.G3VectorQuat([q]), core.G3VectorQuat([p])))
alpha, delta = list(zip(angles[0], angles[1]))
assert(np.allclose(maps.ang_to_quat(alpha, delta), core.G3VectorQuat([maps.c_ang_to_quat_(a[0], a[1]) for a in angles])))

assert(np.allclose(maps.quat_to_ang(maps.c_ang_to_quat_(angle[0], angle[1])), angle))

start = core.G3Time(0)
stop = core.G3Time(10)
q_timestream = maps.ang_to_quat(alpha, delta, start, stop)
assert q_timestream.start == start
assert q_timestream.stop == stop
