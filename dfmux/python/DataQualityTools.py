from spt3g import core
from spt3g.dfmux import DfMuxMetaSample
import numpy as np

def get_empty_timepoint(sample_time):
    tp_frame = core.G3Frame(core.G3FrameType.Timepoint)
    tp_frame["EventHeader"] = sample_time
    tp_frame["DfMux"] = DfMuxMetaSample()
    return tp_frame

class AddMissingTimepoints(object):
    def __init__(self, expected_sample_rate = (20e6 * core.G3Units.Hz)/ 2**17):
        self.samp_spacing_ = 1.0 / expected_sample_rate
        self.prev_sample_time_ = None

    def __call__(self, frame):
        if frame.type != core.G3FrameType.Timepoint:
            return            
        if 'EventHeader' not in frame:
            return
        new_time = frame['EventHeader']
        if self.prev_sample_time_ == None:
            self.prev_sample_time_ = new_time
            return
        t_delt = float(new_time - self.prev_sample_time_)
        n_missing_samples = int(round(t_delt / samp_spacing_ )) - 1
        samp_times = [new_time - samp_spacing_ * i for i in range(1, 1+n_missing_samples)]
        new_packets = [get_empty_timepoint(st) for st in samp_times]
        new_packets.append(frame)
        return new_packets

class FillMissingTimepointFrames(object):
    '''
    Search for Timepoint frames that come at larger-than-expected intervals.
    Fill in the missin frames with garbage, and mark them as garbage for
    conversion to NaNs later.
    '''
    def __init__(self):
        self.dt = None
        self.last_frame = None
        self.buffer = []

    def __call__(self, fr):
        if fr.type != core.G3FrameType.Timepoint:
            return
        if self.dt is None:
            # Calculate the sample rate from the first 20 frames
            self.buffer.append(fr)
            if len(self.buffer) < 20:
                return []
            else:
                times = [_fr['EventHeader'].time for _fr in self.buffer]
                intervals = np.diff(times)
                self.dt = min(intervals[intervals > 0])
                # Now process the frames in the buffer
                out = sum([self._find_and_fill_frames(_fr) for _fr in self.buffer], [])
                self.buffer = []
                return out
        return self._find_and_fill_frames(fr)

    def _find_and_fill_frames(self, fr):
        # This function does most of the work:
        # Find frames with separation that is too large, and fill the gaps
        # with new frames.
        if self.last_frame is None:
            self.last_frame = fr
            return [fr]
        this_dt = fr['EventHeader'].time - self.last_frame['EventHeader'].time
        if this_dt > 1.5 * self.dt:
            # We're missing at least one frame
            nframes = int(np.round(this_dt / self.dt) - 1)
            injected_dt = this_dt / (nframes + 1)
            injected_frs = [core.G3Frame(core.G3FrameType.Timepoint)
                            for i in range(nframes)]
            for i, new_fr in enumerate(injected_frs):
                new_fr['EventHeader'] = self.last_frame['EventHeader'] + \
                                     injected_dt * (i + 1)
                new_fr['DfMux'] = fr['DfMux'].copy()
                try:
                    new_fr['CalibratorOn'] = fr['CalibratorOn']
                except KeyError:
                    pass
                new_fr['GarbageData'] = True
            self.last_frame = fr
            return injected_frs + [fr]
        self.last_frame = fr
        return [fr]

