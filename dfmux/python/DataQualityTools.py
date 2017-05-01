from spt3g import core
from spt3g.dfmux import DfMuxMetaSample

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

