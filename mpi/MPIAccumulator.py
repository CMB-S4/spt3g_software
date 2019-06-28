from spt3g import core

class MPIAccumulator(object):
    '''
    Accumulate data from many frames in memory, sharing metadata for all the
    frames between nodes in an MPI communicator, potentially organized into
    groups (for example by observation ID).

    The result will be stored in the member variable 'fullobs' after
    pipeline termination. 'fullobs' is a dictionary containing entries
    for each data group (see the extractfunc parameter), which in turn
    are lists of 3-element tuples with the following content:
    - First element is the metadata for the frame (see extractfunc)
    - Second element is the node on which the full data exist
    - Third element is the data from the frame (see extractfunc)
      or None if the data are on a remote process.
    '''
    def __init__(self, mpicomm, extractfunc=None, sorter=None, dataframes=[core.G3FrameType.Scan, core.G3FrameType.Timepoint]):
        '''
        Data will be shared between nodes participating in the communicator
        mpicomm. The data stored are based on frames of the types in the
        dataframes argument and are organized based on the results of
        extractfunc and sorter.

        extractfunc is expected to return a three-element tuple:
        - First element is a group ID for the data (for example, an observation
          ID). This is used only for organization and can be anything that
          Python can compare. By default, this is SourceName-ObservationID.
        - Second element is metadata. This information is shared between
          all processes participating in the communicator at completion and
          should include any information you need to assess whether the data
          are interesting for future calculations (for example, the start time
          of the frame). By default, frames are sorted by this value. None
          by default.
        - Third element is the data. This stays node-local. It could be just
          a copy of the frame (the default) or just the information (in any
          format understandable to Python) from the frame you expect to use
          later.
        Note that extractfunc is called on all frame types (including
        calibration data) but *only the results from frames in dataframes*
        are stored. Thus, it is the responsibility of extractfunc to cache
        any applicable calibration, etc. information and store it in the
        return value for data frames if needed.

        sorter is used, if defined, to sort the list of frame data in each
        group. This is passed to the Python sorted() function as a 'key'
        argument and receives the three-element tuple stored for each set of
        frame data (metadata, node, data). By default, does no sorting
        (frames will still appear in a common order across nodes related to
        the data distribution if unsorted).
        '''
        self.mpicomm = mpicomm
        self.extractfunc = extractfunc
        def defaultextract(f):
            if f.type not in self.dataframes:
                return
            return ('%s-%d' % (f['SourceName'], f['ObservationID']), None, f)
        if self.extractfunc is None:
            self.extractfunc = defaultextract
        self.sorter = sorter
        self.dataframes = dataframes
        self.localdata = {}
    def __call__(self, frame):
        if frame.type == core.G3FrameType.EndProcessing:
            self.finalize()
            return

        if frame.type not in self.dataframes:
            self.extractfunc(frame) # Give it a chance to cache things
            return

        obskey, metadata, data = self.extractfunc(frame)
        if obskey not in self.localdata:
            self.localdata[obskey] = []
        self.localdata[obskey].append((metadata, data))

    def finalize(self):
        # Put all observation metadata everywhere
        fullobs = {}
        localmetadata = {}
        for k,v in self.localdata.items():
            localmetadata[k] = [(m, self.mpicomm.rank, None) for m,d in v]
        for i in range(self.mpicomm.size):
            if i == self.mpicomm.rank:
                self.mpicomm.bcast(localmetadata, root=i)
                metadata = {}
                for k,v in self.localdata.items():
                    metadata[k] = [(m, i, d) for m,d in v]
            else:
                metadata = self.mpicomm.bcast(None, root=i)
            for obs, data in metadata.items():
                if obs not in fullobs:
                    fullobs[obs] = []
                fullobs[obs] += data
        del localmetadata
        del self.localdata
            
        # Now let's organize it
        if self.sorter is not None:
            for obs, data in fullobs.items():
                data.sort(key=self.sorter)

        self.fullobs = fullobs

