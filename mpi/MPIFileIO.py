import mpi4py
from spt3g import core

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
    will be read and processed by the remaining nodes.
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
 
