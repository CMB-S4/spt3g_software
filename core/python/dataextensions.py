from spt3g.core import G3Bool, G3Int, G3Double, G3String

__all__ = []

def g3data_comparison_wrap(cls):
    cls.__eq__ = lambda a, b: a.value == b.value
    cls.__neq__ = lambda a, b: a.value != b.value
    cls.__ge__ = lambda a, b: a.value >= b.value
    cls.__gt__ = lambda a, b: a.value > b.value
    cls.__le__ = lambda a, b: a.value <= b.value
    cls.__lt__ = lambda a, b: a.value < b.value

g3data_comparison_wrap(G3Bool)
g3data_comparison_wrap(G3Int)
g3data_comparison_wrap(G3Double)
g3data_comparison_wrap(G3String)
