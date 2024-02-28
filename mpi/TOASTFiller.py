from spt3g import core, calibration, dfmux
import numpy
from .MPIAccumulator import MPIAccumulator
import pickle

class TOASTFiller(object):
    '''
    Fill a TOAST TOD (stored in the toastobses member after completion) from
    a frame stream. Data will be distributed across the given MPI communicator
    in the way that TOAST would like.
    '''
    def __init__(self, mpicomm, boresight='OnlineRaDecRotation', timestreams='RawTimestreams_I', wiring='WiringMap'):
        self.boresight = boresight
        self.timestreams = timestreams
        self.wiring = wiring
        self.mpicomm = mpicomm
        self.mpia = MPIAccumulator(mpicomm, self.extractor, lambda m: m[0]['start'], dataframes=[core.G3FrameType.Scan])

    def extractor(self, frame):
        if frame.type == core.G3FrameType.Wiring:
            self.keys = frame[self.wiring].keys()
            self.invkeys = {k: i for i,k in enumerate(self.keys)}

        if frame.type == core.G3FrameType.Calibration:
            # from toast import qarray as qa
            #self.detquat = {k: qa.from_angles(v.y_offset/core.G3Units.rad, v.x_offset/core.G3Units.rad, 0) for k,v in frame['BolometerProperties'].items()}
            self.detquat = None

        if frame.type == core.G3FrameType.Scan:
            timestreams = frame[self.timestreams]
            bs = numpy.asarray(frame[self.boresight])
            bs = bs[:,[1,2,3,0]] # Swizzle for TOAST
            scan = {'dets': self.keys, 'boresight': bs, 'detquat': self.detquat, 'obsid': '%s-%d' % (frame['SourceName'], frame['ObservationID']), 'start': timestreams.start, 'stop': timestreams.stop, 'nsamples': len(timestreams.values()[0])}
            return scan['obsid'], scan, timestreams


    def __call__(self, frame):
        rv = self.mpia(frame)
        if frame.type != core.G3FrameType.EndProcessing:
            return rv

        from toast.tod import TODCache

        # Build TOAST TOD from our full data set, now that we have it
        toastobses = []
        for obs,scans in self.mpia.fullobs.items():
            nsamptot = sum([i[0]['nsamples'] for i in scans])
            d = TODCache(self.mpicomm, scans[0][0]['dets'], nsamptot, detquats=scans[0][0]['detquat'])

            t = numpy.concatenate([numpy.linspace(i[0]['start'].time/core.G3Units.s, i[0]['stop'].time/core.G3Units.s, i[0]['nsamples']) for i in scans])
            d.write_times(local_start=0, stamps=t[slice(*d.local_samples)])
            bs = numpy.concatenate([i[0]['boresight'] for i in scans])
            d.write_boresight(local_start=0, data=bs[slice(*d.local_samples)])

            startsample = 0
            for i in scans:
                i[0]['startsamp'] = startsample
                startsample += i[0]['nsamples']

            # Make a list of who needs what
            dataneeds = (self.mpicomm.rank, d.local_samples, d.local_dets)
            dataneeds = self.mpicomm.allgather(dataneeds)

            # Send/receive needed data
            outdata = [[] for i in range(self.mpicomm.size)]
            for need in dataneeds:
                for i in scans:
                    if i[2] is None:
                        continue
                    mask = slice(max(need[1][0], i[0]['startsamp']),
                            min(need[1][1] + need[1][0], i[0]['startsamp'] +
                                i[0]['nsamples']))
                    if mask.stop <= mask.start:
                        continue # No overlap

                    # Scan overlaps a need, get the relevant timestreams
                    maskstart = mask.start
                    mask = slice(mask.start - i[0]['startsamp'], mask.stop - i[0]['startsamp'])
                    chunk = {}
                    for k in need[2]:
                        chunk[k] = numpy.asarray(i[2][k])[mask]
                    outdata[need[0]].append((maskstart, chunk))
            del dataneeds

            outdata = self.mpicomm.alltoall(outdata) # Swap with everyone

            # Now stitch everything into the TOAST structure
            for sourcenode in outdata:
                for chunk in sourcenode:
                    localstart = chunk[0] - d.local_samples[0]
                    for k,v in chunk[1].items():
                        d.write(k, localstart, v)

            del self.mpia # Don't need this anymore
            toastobses.append({'id': obs, 'tod': d})

        self.toastobses = toastobses
        return frame

