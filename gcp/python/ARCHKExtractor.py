import numpy
from spt3g import core, dfmux

@core.indexmod
class UnpackSPTpolHKData(object):
    '''
    Extracts SPTpol housekeeping information from ARC files
    '''

    def __init__(self):
        self.hk_seen = False

    def __call__(self, f):
        if f.type != core.G3FrameType.GcpSlow:
            return

        hkregs = f['array']['muxHousekeeping']

        # HK registers are just duplicated when unchanged. Deduplicate them here
        if self.hk_seen and not numpy.asarray(hkregs['hkIsValid']).any():
            return
        self.hk_seen = True

        hkout = dfmux.DfMuxHousekeepingMap()

        for i, serial in enumerate(hkregs['moboId']):
            board = dfmux.HkBoardInfo()
            board.timestamp = hkregs['utc'] # TimeStamp is off by exactly 30 days, so ignore it
            board.timestamp_port = 'UNKNOWN'
            board.serial = str(serial)
            board.fir_stage = hkregs['firStage'][i]
            board.voltages['Ext5V5N'] = hkregs['mbExt5V5N'][i]
            board.voltages['Ext5V5P'] = hkregs['mbExt5V5P'][i]
            board.voltages['Int5V5N'] = hkregs['mbInt5V5N'][i]
            board.voltages['Int5V5P'] = hkregs['mbInt5V5P'][i]
            board.voltages['Squid8V5N'] = hkregs['mbSquid8V5N'][i]
            board.voltages['Squid8V5P'] = hkregs['mbSquid8V5P'][i]
            board.voltages['VAuxN'] = hkregs['mbVAuxN'][i]
            board.voltages['VAuxP'] = hkregs['mbVAuxP'][i]
            # XXX: mbSensorX?

            for j, mezz in enumerate(['A', 'B']):
                mezzhk = dfmux.HkMezzanineInfo()
                mezzhk.power = True
                mezzhk.present = True
                mezzhk.serial = 'N/A'
                mezzhk.part_number = 'SPTpol Mezz'
                mezzhk.revision = ''
    
                for v in ['V3P', 'V5N', '_VDemodA', '_VDemodB']:
                    mezzhk.voltages[v.strip('_')] = hkregs['mezz' + mezz + v][i]

                for k, mod in enumerate(['A', 'B']):
                    modhk = dfmux.HkModuleInfo()
                    modhk.module_number = k
                    modhk.carrier_gain = hkregs['carrierGain'][i][j*2 + k]
                    modhk.nuller_gain = hkregs['nullerGain'][i][j*2 + k]
                    modhk.demod_gain = hkregs['demodGain'][i][j*2 + k]
                    if 'mezz' + mezz + 'Car' + mod + 'Overload' in hkregs: 
                        modhk.carrier_railed = hkregs['mezz' + mezz + 'Car' + mod + 'Overload'][i]
                        modhk.nuller_railed = hkregs['mezz' + mezz + 'Nul' + mod + 'Overload'][i]
                    modhk.demod_railed = hkregs['mezz' + mezz + 'Adc' + mod + 'Overload'][i]

                    nchan = len(hkregs['carrierSettings'][0])//(3*2*2) # 3 fields per mezz/module
                    for chan in range(nchan):
                        chanhk = dfmux.HkChannelInfo()
                        chanhk.channel_number = chan+1
                        arcchan = (j*2 + k)*nchan + chan

                        if 'danGain' in hkregs:
                            chanhk.dan_gain = hkregs['danGain'][i][arcchan]
                            chanhk.dan_feedback_enable = hkregs['danFeedbackEnable'][i][arcchan]
                            chanhk.dan_streaming_enable = hkregs['danStreamingEnable'][i][arcchan]
                            chanhk.dan_accumulator_enable = hkregs['danAccumulatorRailed'][i][arcchan]
                            chanhk.dan_railed = hkregs['danAccumulatorRailed'][i][arcchan]

                        # Decrypt xSettings registers: (freq, phase, amp)
                        chanhk.carrier_frequency = hkregs['carrierSettings'][i][arcchan*3 + 0]*25.e6/(2.**32-1)*core.G3Units.Hz
                        chanhk.carrier_amplitude = hkregs['carrierSettings'][i][arcchan*3 + 2]/(2.**23 - 1.)
                        chanhk.nuller_amplitude = hkregs['nullerSettings'][i][arcchan*3 + 2]/(2.**23 - 1.)
                        chanhk.demod_frequency = hkregs['demodSettings'][i][arcchan*2 + 0]*25.e6/(2.**32-1)*core.G3Units.Hz

                        modhk.channels[chan+1] = chanhk

                    mezzhk.modules[k+1] = modhk

                board.mezz[j+1] = mezzhk

            hkout[serial] = board

        hkframe = core.G3Frame(core.G3FrameType.Housekeeping)
        hkframe['DfMuxHousekeeping'] = hkout

        return [f, hkframe]

