#!/usr/bin/env python
from spt3g import core, maps
import numpy as np

a = core.quat(2,3,4,5)
a2 = core.quat(2,1,8,5)
b = core.G3VectorQuat([a, a**2, 2*a, a2])

assert(np.isclose(maps.quat_to_ang(a), maps.c_quat_to_ang_(a)).all())
assert(np.isclose(maps.quat_to_ang(b), np.asarray([maps.c_quat_to_ang_(x) for x in b]).transpose()).all())

angle = (.3, 0.4)
angles = ((.3, 0.8), (0.4, 0.2))

assert(maps.ang_to_quat(angle[0], angle[1]) == maps.c_ang_to_quat_(angle[0], angle[1]))
alpha, delta = list(zip(angles[0], angles[1]))
assert((np.asarray(maps.ang_to_quat(alpha, delta)) == np.asarray(core.G3VectorQuat([maps.c_ang_to_quat_(a[0], a[1]) for a in angles]))).all())

assert(np.isclose(maps.quat_to_ang(maps.c_ang_to_quat_(angle[0], angle[1])), angle).all())

