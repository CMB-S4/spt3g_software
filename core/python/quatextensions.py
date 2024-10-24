import numpy as np
from . import Quat, G3VectorQuat, G3TimestreamQuat

__all__ = []

quat_types = (Quat, G3VectorQuat, G3TimestreamQuat)


def quat_ufunc(self, ufunc, method, *inputs, **kwargs):
    """Numpy ufunc interface for vectorized quaternion methods."""
    if ufunc.__name__ in ["isinf", "isnan", "isfinite"] and len(inputs) == 1:
        v = getattr(ufunc, method)(np.asarray(inputs[0]), **kwargs)
        if ufunc.__name__ == "isfinite":
            return np.all(v, axis=-1)
        return np.any(v, axis=-1)
    if ufunc.__name__.startswith("logical"):
        args = []
        for arg in inputs:
            if isinstance(arg, quat_types):
                arg = np.any(np.asarray(arg), axis=-1)
            args.append(arg)
        return getattr(ufunc, method)(*args, **kwargs)
    if method != "__call__" or kwargs:
        return NotImplemented
    if len(inputs) == 1:
        if ufunc.__name__ == "absolute":
            return self.abs()
        if ufunc.__name__ == "negative":
            return self.__neg__()
        if ufunc.__name__ == "conjugate":
            return self.conj()
        if ufunc.__name__ == "reciprocal":
            return Quat(1, 0, 0, 0) / self
    if len(inputs) == 2 and np.isscalar(inputs[1]):
        if ufunc.__name__ == "power":
            return self.__pow__(inputs[1])
    return NotImplemented


Quat.__array_ufunc__ = quat_ufunc
G3VectorQuat.__array_ufunc__ = quat_ufunc
G3TimestreamQuat.__array_ufunc__ = quat_ufunc
