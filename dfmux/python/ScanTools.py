from spt3g import core

@core.indexmod
class FixedLengthScans(object):
    '''Makes scans of length N timepoints.'''
    def __init__(self, N=1000):
        self.N = N
        self.count = 0
    def __call__(self, frame):
        ret = []
        if frame.type == core.G3FrameType.Timepoint:
            if self.count % self.N == 0:
                ret.append(core.G3Frame(core.G3FrameType.Scan))
            self.count += 1
        ret.append(frame)
        return ret

