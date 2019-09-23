#!/usr/bin/env python
from spt3g import core, coordinateutils
import numpy as np

a = core.quat(2,3,4,5)
a2 = core.quat(2,1,8,5)
b = core.G3VectorQuat([a, a**2, 2*a, a2])

assert(np.isclose(coordinateutils.quat_to_ang(a), coordinateutils.c_quat_to_ang(a)).all())
assert(np.isclose(coordinateutils.quat_to_ang(b), np.asarray([coordinateutils.c_quat_to_ang(x) for x in b]).transpose()).all())

angle = (.3, 0.4)
angles = ((.3, 0.8), (0.4, 0.2))

assert(coordinateutils.ang_to_quat(angle[0], angle[1]) == coordinateutils.c_ang_to_quat(angle[0], angle[1]))
alpha, delta = list(zip(angles[0], angles[1]))
assert((np.asarray(coordinateutils.ang_to_quat(alpha, delta)) == np.asarray(core.G3VectorQuat([coordinateutils.c_ang_to_quat(a[0], a[1]) for a in angles]))).all())

assert(np.isclose(coordinateutils.quat_to_ang(coordinateutils.c_ang_to_quat(angle[0], angle[1])), angle).all())

