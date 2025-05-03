from . import G3Frame, G3FrameType

__all__ = []

@staticmethod
def str_to_frame_types(types):
    """
    Takes a string of characters representing frame types,
    and returns a list of the corresponding G3FrameType subclasses.
    """
    frame_types = []
    for typenum in map(ord,types):
        if (typenum in G3FrameType.values):
            frame_types.append(G3FrameType.values.get(typenum))
        else:
            raise ValueError("Invalid Frame Type: {0}".format(
                             chr(typenum)))
    return frame_types

def items(self):
    """
    items implementation for G3Frame
    """
    keys = self.keys()
    for key in keys:
        yield (key, self[key])

def __iter__(self):
    """
    __iter__ implementation for G3Frame
    """
    keys = self.keys()
    for key in keys:
        yield key

G3FrameType.from_string = str_to_frame_types
G3Frame.items = items
G3Frame.__iter__ = __iter__
