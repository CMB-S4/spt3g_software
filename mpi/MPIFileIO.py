from mpi4py import MPI
from spt3g import core
import random

class MPIFileReader(object):
    '''
    Do parallel I/O and processing across an MPI communicator. The style of
    parallelism here is that each process in the communicator gets a file, reads
    the file, and processes the file. Supports shared files, which are read
    from the lead process and then broadcasted to the others before they read
    their own.

    The list of files is specified as a list of lists. For example,
    ['offline_calibration.g3', ['0000.g3', '0002.g3', '0003.g3']]
    will cause all nodes to process 'offline_calibration.g3' and then the
    remaining files (in the parallel block denoted by the nested list)
    will be read and processed by the remaining nodes in a distributed way.
    '''
    def __init__(self, mpicomm, files):
        self.comm = mpicomm
        self.files = files
        self.cur_reader = None
        self.broadcast = False
    def nextreader(self):
        nextfilelist = self.files.pop(0)
        self.broadcast = False
        if isinstance(nextfilelist, str):
            if self.comm.rank == 0:
                self.broadcast = True
                self.cur_reader = core.G3Reader(nextfilelist)
            else:
                self.cur_reader = lambda fr: self.comm.bcast(None,root=0)
        else:
            toreadhere = nextfilelist[self.comm.rank::self.comm.size]
            self.cur_reader = core.G3Reader(toreadhere)
    def __call__(self, frame):
        assert(frame is None)
        out = []
        while len(out) == 0:
            if self.cur_reader is None:
                if len(self.files) == 0:
                    return []
                self.nextreader()
            out = self.cur_reader(frame)
            if self.broadcast:
                self.comm.bcast(out,root=0)
            if len(out) == 0:
                # Move to the next reader if this one is done
                self.cur_reader = None
        return out

class MPIFrameParallelizer(object):
    '''
    Do parallel I/O and processing across an MPI communicator. The style of
    parallelism here is that a set of IO processes read files from disk
    and distribute frames round-robin-style across a worker communicator.

    Rules to make this work:
    - First module on IO nodes must be MPIFileReader
    - Second module on IO nodes must be DeduplicateMetadata
    - Third module on IO nodes, first on workers, must be MPIFrameParallelizer

    You may want to use the MPIIODistributor pipe segment to make sure these
    rules are followed.
    '''
    def __init__(self, iocomm, cpucomm, cpucomm_startrank, dataframetype=[core.G3FrameType.Timepoint, core.G3FrameType.Scan]):
        self.iocomm = iocomm
        self.cpucomm = cpucomm
        self.cpucomm_startrank = cpucomm_startrank
        self.dataframes = dataframetype
    def __call__(self, frame):
        if frame is None:
            # Initial frame is None on CPU nodes
            frame = self.cpucomm.recv()
            if frame.type == core.G3FrameType.EndProcessing:
                frame = [] # This is the start module
            return frame
        if frame.type not in self.dataframes:
            # Synchronize and broadcast metadata; should be the same
            # in all streams
            self.iocomm.barrier()
            if self.iocomm.rank == 0:
                checkframe = self.iocomm.bcast(frame, root=0)
            else:
                checkframe = self.iocomm.bcast(None, root=0)
            assert(checkframe.__getstate__()[0] == frame.__getstate__()[0]) # All metadata the same?
            # Get metadata frame in the queue on every worker process
            for i in list(range(self.cpucomm_startrank, self.cpucomm.size))[self.iocomm.rank::self.iocomm.size]:
                self.cpucomm.send(frame, i)
        else:
                # Send it somewhere random
                self.cpucomm.send(frame, random.randint(self.cpucomm_startrank, self.cpucomm.size-1))

        if frame.type == core.G3FrameType.EndProcessing:
            return frame

        return [] # Terminate processing on IO nodes
 
@core.pipesegment
def MPIIODistributor(pipe, mpicomm=MPI.COMM_WORLD, n_io=10, files=[]):
    subcomm = mpicomm.Split(mpicomm.rank < n_io, mpicomm.rank)
    if mpicomm.rank < n_io:
        pipe.Add(MPIFileReader, mpicomm=subcomm, files=files)
        pipe.Add(core.util.DeduplicateMetadata)
    pipe.Add(MPIFrameParallelizer, iocomm=(subcomm if mpicomm.rank < n_io else None), cpucomm=mpicomm, cpucomm_startrank=n_io)

