from spt3g.core import G3Frame, G3FrameType

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

@property
def type_enum_key(self):
    """
    Return G3FrameType's C++ enum character key.
    """
    return chr(self.real)

def iteritems(self):
    """
    iteritems implementation for G3Frame
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

def get(self, key, default=None):
    """
    get implementation for G3Frame
    """
    if key in self.keys():
        return self[key]
    return default

def setdefault(self, key, default=None):
    """
    setdefault implementation for G3Frame
    """
    if key not in self.keys():
        self[key] = default

def pop(self, key, default=None):
    """
    pop implementation for G3Frame
    """
    if key not in self.keys():
        return default
    value = self[key]
    del self[key]
    return value

G3FrameType.from_string = str_to_frame_types
G3FrameType.key = type_enum_key
G3Frame.iteritems = iteritems
G3Frame.items = iteritems
G3Frame.__iter__ = __iter__
G3Frame.get = get
G3Frame.setdefault = setdefault
G3Frame.pop = pop

del str_to_frame_types
del type_enum_key
del iteritems
del __iter__
del get
del setdefault
del pop
