import numpy as np
from . import Quat, G3VectorQuat, G3TimestreamQuat

__all__ = []


def quat_ufunc(self, ufunc, method, *inputs, **kwargs):
    """Numpy ufunc interface for vectorized quaternion methods."""
    if method == "__call__" and len(inputs) == 1 and not kwargs:
        if ufunc.__name__ == "absolute":
            return self.abs()
        if ufunc.__name__ == "conjugate":
            return self.conj()
        if ufunc.__name__ == "reciprocal":
            return Quat(1, 0, 0, 0) / self
    return NotImplemented


Quat.__array_ufunc__ = quat_ufunc
G3VectorQuat.__array_ufunc__ = quat_ufunc
G3TimestreamQuat.__array_ufunc__ = quat_ufunc
