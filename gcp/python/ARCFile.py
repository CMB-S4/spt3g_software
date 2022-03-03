from spt3g.gcp import ARCFileReader, ARCExtract

class ARCFile(object):
    '''Iterable class for ARC files, as created by GCP. Loop through frames by doing something like:
    f = gcp.ARCFile('/path/to/arc.dat')
    for frame in f:
        print( frame )

    An entire file can also be read into an indexable list by doing:
    f = list(gcp.ARCFile('/path/to/arc.dat'))
    '''
    def __init__(self, path, extract=False):
        self.reader = ARCFileReader(path)
        self.extract = extract

    def __iter__(self):
        return self

    def next(self):
        frames = self.reader.Process(None)
        if len(frames) == 0:
            raise StopIteration('No more frames in file')
        if len(frames) > 1:
            raise ValueError('Too many frames returned by reader')
        frame = frames[0]

        if self.extract:
            # calibrate and parse arc frames
            from spt3g.core import G3Pipeline, G3InfiniteSource, InjectFrame

            pipe = G3Pipeline()
            pipe.Add(G3InfiniteSource, n=1)
            pipe.Add(InjectFrame, frame=frame)
            pipe.Add(ARCExtract)
            pipe.Run()

        return frame

    __next__ = next
