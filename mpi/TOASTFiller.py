from toast.tod import TODCache
from toast import qarray as qa
from spt3g import core, calibration, dfmux, coordinateutils
import numpy
from .MPIAccumulator import MPIAccumulator

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
            outdata = [[]]*self.mpicomm.size

            # Send/receive needed data
            for need in dataneeds:
                for i in scans:
                    if i[2] is None:
                        continue
                    if (i[0]['startsamp'] + i[0]['nsamples']) < need[1][0]:
                        continue
                    if i[0]['startsamp'] > need[1][1]:
                        break

                    # Scan overlaps a need, get the relevant timestreams
                    mask = slice(need[1][0] - i[0]['startsamp'], need[1][1] + 1 - i[0]['startsamp'])
                    if mask.start < 0:
                        maskstart = need[1][0]
                        mask = slice(0, mask.stop)
                    else:
                        maskstart = i[0]['startsamp']
                    if mask.stop > i[0]['nsamples']:
                        mask = slice(mask.start, i[0]['nsamples'])
                    chunk = core.G3TimestreamMap()
                    for k in need[2]:
                        chunk[k] = i[2][k][mask]
                    outdata[need[0]].append((maskstart, chunk))
            del self.mpia # Avoid RAM balloon by deleting the old copies
            del dataneeds

            outdata = self.mpicomm.alltoall(outdata) # Swap with everyone

            # Now stitch everything into the TOAST structure
            for chunk in outdata:
                localstart = chunk[0] - d.local_samples[0]
                for k,v in chunk[1].items():
                    d.write(k, localstart, v)

            toastobses.append({'id': obs[0]['obsid'], 'tod': d})

        self.toastobses = toastobses
        return frame

