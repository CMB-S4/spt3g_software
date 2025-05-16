from . import G3FrameType

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

G3FrameType.from_string = str_to_frame_types
