from spt3g.core import indexmod, pipesegment, G3FrameType, log_fatal, G3Reader, G3NetworkSender

@indexmod
def Delete(frame, keys=[], type=None):
    '''
    Deletes specified keys from frame. If optional type specified, only acts on frames of the given type.
    '''
    if type is not None and frame.type != type:
        return

    for key in keys:
        if key in frame:
            del frame[key]

@indexmod
def Rename(frame, keys={}, type=None):
    '''
    Renames specified keys in frame. If optional type specified, only acts on frames of the given type. Argument is a dictionary mapping old names to new ones.
    '''
    if type is not None and frame.type != type:
        return

    for key in keys:
        if key in frame:
            frame[keys[key]] = frame[key]
            del frame[key]

@indexmod
def Dump(frame, type=None, added_message = None):
    '''
    Prints frames to console. If optional type specified, only acts on frames of the given type.
    '''
    if type is not None and frame.type != type:
        return
    if frame.type == G3FrameType.EndProcessing:
        return
    if added_message:
        print(added_message)
    print(frame)


@indexmod
class InjectFrame(object):
    """
    Inject an arbitrary frame into a pipeline.

    Arguments
    ---------
    frame : G3Frame
        The frame to inject
    """
    def __init__(self, frame):
        self.frame = frame
    def __call__(self, frame):
        if self.frame is None:
            return
        out = [self.frame, frame]
        self.frame = None
        return out


@indexmod
def InjectDebug(frame, type=None, debug_start_func = None):
    '''starts a pdb session when a frame of type shows up.

    The frame data is stored in the variable names "frame".

    If debug_start_func is not None, only starts a debug session when 
            debug_start_func(frame) == True
        '''
    if type is None or frame.type == type:
        if ((debug_start_func is None) or debug_start_func(frame)):
            import pdb, rlcompleter
            pdb.Pdb.complete = rlcompleter.Completer(locals()).complete
            pdb.set_trace()


@indexmod
class AbortAfterNFrames(object):
    '''Stops processing after n_frames frames go by'''
    def __init__(self, type, n_frames):
        self.n_desired_frames = n_frames
        self.num_frames = 0
        self.type = type
    def __call__(self, frame):
        if self.num_frames >= self.n_desired_frames:
            log_fatal("Manual Abort Triggered")
        if frame.type == self.type:
            self.num_frames += 1
        return

@pipesegment
def G3NetworkReceiver(pipe, hostname='localhost', port=5978):
    '''
    Emulation of old G3NetworkReceiver class. Equivalent to pointing
    G3Reader at a TCP URL.
    '''

    pipe.Add(G3Reader, filename='tcp://' + hostname + ':' + str(port))

@indexmod
class G3ThrottledNetworkSender(object):
    '''
    Send every Nth frame of certain types using a wrapped G3NetworkSender.
    All instances of frames not in the dictionary frame_decimation will be sent
    at their full rate.
    '''
    def __init__(self, hostname='*', port=5978, frame_decimation = {G3FrameType.Timepoint: 10}, max_queue_size=0):
        self.sender = G3NetworkSender(hostname=hostname, port=port, max_queue_size=max_queue_size)
        self.decimation = frame_decimation
        self.counts = {}
        for k in self.decimation.keys():
            self.counts[k] = 0
    def __call__(self, frame):
        if frame.type in self.counts:
            if self.decimation[frame.type] == 0:
                return
            self.counts[frame.type] += 1
            if self.counts[frame.type] % self.decimation[frame.type] != 0:
                return

        self.sender(frame)

@indexmod
class DeduplicateMetadata(object):
    '''
    Drop metadata frames (e.g. Calibration, Wiring) for which the previous
    metadata frame of the same type is byte-for-byte identical. This can be
    handy when, for example, reading in many files from the G3MultiFileWriter,
    which copies metadata frames to the beginning of each file. Considers
    all frames not in <dataframetypes> to be metadata (by default, everything
    except Timepoint and Scan frames).
    '''
    def __init__(self, dataframetype=[G3FrameType.Timepoint, G3FrameType.Scan]):
        self.dataframes = dataframetype
        self.metacache = {}
    def __call__(self, f):
        if f.type in self.dataframes or f.type == G3FrameType.EndProcessing:
            # Pass data or end processing frames through
            return

        # Compare serialized form of frame
        from hashlib import md5
        rawdata = md5(f.__getstate__()[1]).hexdigest()
        if f.type in self.metacache and self.metacache[f.type] == rawdata:
            # Same as the last one, so drop this frame
            return []

        # Otherwise, pass through and cache
        self.metacache[f.type] = rawdata

@indexmod
class DropOrphanMetadata(object):
    '''
    Remove metadata frames (e.g. Calibration, Wiring) without
    intervening data frames (e.g. Timepoint, Scan, specified by the
    <dataframetype> argument to the constructor). The metadata frames that
    do show up will be the most recent of each type and appear in their
    original order.
    '''
    def __init__(self, dataframetype=[G3FrameType.Timepoint, G3FrameType.Scan]):
        self.dataframes = dataframetype
        self.metacache = []
    def __call__(self, frame):
        if frame.type == G3FrameType.EndProcessing:
            return
        if frame.type in self.dataframes:
            # If we got a data frame, dump any queued metadata before it
            # and clear the queue.
            out = self.metacache + [frame]
            self.metacache = []
            return out
        else:
            # Insert new metadata frames into the same spot in the queue
            # as the last instance of a frame of that type or the end of the
            # queue if there is no such frame.
            for i in range(len(self.metacache)):
                if self.metacache[i].type == frame.type:
                    self.metacache[i] = frame
                    break
            else:
                self.metacache.append(frame)
            return []
                    

del indexmod
del pipesegment
