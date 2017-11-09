from spt3g.core import G3FrameObject

__all__ = []

def get_value(obj):
    if isinstance(obj, G3FrameObject):
        return obj.value
    return obj

def g3data_comparison_wrap(cls):
    cls.__eq__ = lambda a, b: a.value == get_value(b)
    cls.__neq__ = lambda a, b: a.value != get_value(b)
    cls.__ge__ = lambda a, b: a.value >= get_value(b)
    cls.__gt__ = lambda a, b: a.value > get_value(b)
    cls.__le__ = lambda a, b: a.value <= get_value(b)
    cls.__lt__ = lambda a, b: a.value < get_value(b)

g3data_comparison_wrap(G3FrameObject)
