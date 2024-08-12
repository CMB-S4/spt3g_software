from spt3g import core

__all__ = ["G3File"]


class G3File(object):
    """Iterable class for G3 files, as created by G3Writer. Loop through frames by doing something like:

        with core.G3File("/path/to/file.g3") as f:
            for frame in f:
                print(frame)

    An entire file can also be read into an indexable list by doing:

        f = list(core.G3File("/path/to/file.g3"))
    """

    def __init__(self, path):
        self.reader = core.G3Reader(path)

    def __iter__(self):
        return self

    def next(self):
        if self.reader is None:
            raise StopIteration("Reader closed")
        frames = self.reader.Process(None)
        if len(frames) == 0:
            raise StopIteration("No more frames in file")
        if len(frames) > 1:
            raise ValueError("Too many frames returned by reader")
        return frames[0]

    __next__ = next

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        del self.reader
        self.reader = None


# writer context methods


def writer_enter(self):
    return self


def writer_exit(self, *args, **kwargs):
    fr = core.G3Frame(core.G3FrameType.EndProcessing)
    self(fr)


core.G3Writer.__enter__ = writer_enter
core.G3Writer.__exit__ = writer_exit
core.G3MultiFileWriter.__enter__ = writer_enter
core.G3MultiFileWriter.__exit__ = writer_exit
