from spt3g import core
from spt3g.dfmux import DfMuxHousekeepingMap, HkBoardInfo, HkMezzanineInfo, HkModuleInfo, HkChannelInfo
from .TuberClient import TuberClient
import socket, struct, time

@core.indexmod
class LegacyHousekeepingConsumer(object):
    '''
    Collect housekeeping data from the legacy (i.e. SPTpol-ish) mux boards
    defined in the wiring map. Will add a key called 'DfMuxHousekeeping' to
    any housekeeping frame that goes by containing the data as of the arrival
    of the housekeeping frame. Use in conjunction with a
    dfmux.PeriodicHousekeepingCollector to get data at fixed intervals.
   
    If collecting real-time data, you may want to set subprocess=True when
    adding this module.
    '''
    def __init__(self):
        self.tuber = None
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Wiring:
            wmap = frame['WiringMap']
            board_ips = set([m.board_ip for m in wmap.values()])

            self.tuber = [(ip, TuberClient(socket.inet_ntoa(struct.pack('i', m.board_ip)))) for ip in board_ips]

        if frame.type == core.G3FrameType.Housekeeping:
            for board,board_tuber in self.tuber:
                board_tuber.CallMethod('Dfmux', ['classic_dfmux_dump', 'classic_housekeeping'], False)

            hkdata = DfMuxHousekeepingMap()
            for board,board_tuber in self.tuber:
                dat = dict(board_tuber.GetReply()[0]['result'])
                dat.update(board_tuber.GetReply()[1]['result'])

                serial = struct.unpack(">i", struct.pack("i", board))[0] & 0xff
                hkdata[serial] = self.HousekeepingFromJSON(dat)
            frame['DfMuxHousekeeping'] = hkdata

    @classmethod
    def HousekeepingFromJSON(cls, dat):
        '''
        Build HKBoardInfo object from the union of the JSON blobs
        returned by the classic_housekeeping and class_dfmux_dump
        Tuber calls
        '''
        # Board-global quantities
        boardhk = HkBoardInfo()

        if dat['ts']['port'] == 'IRIG test':
            boardhk.timestamp = core.G3Time(dat['ts']['s']*core.G3Units.s)
        else:
            year = dat['ts']['y']
            if year == 0:
                # It probably isn't 1900
                systime = time.gmtime()
                # Check for New Year's, assuming no more than 24 hours clock slew
                if dat['ts']['d'] == 1 and systime.tm_yday >= 365:
                    year = systime.tm_year + 1
                elif dat['ts']['d'] >= 365 and systime.tm_yday == 1:
                    year = systime.tm_year - 1
                else:
                    year = systime.tm_year
                year -= 2000
            boardhk.timestamp = core.G3Time(y=year,d=dat['ts']['d'],h=dat['ts']['h'],m=dat['ts']['m'],s=dat['ts']['s'],ss=dat['ts']['ss'])

        boardhk.timestamp_port = str(dat['ts']['port'])
        boardhk.serial = 'N/A'
        boardhk.fir_stage = dat['fir_stage']
        for i in dat['voltages']['mb'].items():
            boardhk.voltages[str(i[0])] = i[1]
        for i in dat['temperatures']['mb'].items():
            boardhk.temperatures[str(i[0])] = i[1]

        # Mezzanines
        for mezz in [1,2]:
            mezzhk = HkMezzanineInfo()
            mezzhk.power = dat['voltages']['mezz%d' % mezz]['v3p'] > 0
            mezzhk.present = True
            mezzhk.serial = 'N/A'
            mezzhk.part_number = 'SPTpol Mezz'
            mezzhk.revision = 'N/A'
            for i in dat['voltages']['mezz%d' % mezz].items():
                if isinstance(i[1], float):
                    mezzhk.voltages[str(i[0])] = i[1]

            # Modules
            for mod in [1,2]:
                wire = (mezz - 1)*2 + mod

                # XXX: SQUID data (bias point, FLL) not in HK call!
                # XXX: SQUID tuning state
                modhk = HkModuleInfo()
                modhk.carrier_gain = dat['mezz_gains']['wire%d' % wire]['carrier']
                modhk.nuller_gain = dat['mezz_gains']['wire%d' % wire]['nuller']
                modhk.demod_gain = dat['mezz_gains']['wire%d' % wire]['demod']

                modhk.carrier_railed = dat['voltages']['mezz%d' % mezz]['overload_dmfs_car_%s' % (['a', 'b'][mod-1])]
                modhk.nuller_railed = dat['voltages']['mezz%d' % mezz]['overload_dmfs_nul_%s' % (['a', 'b'][mod-1])]
                modhk.demod_railed = dat['voltages']['mezz%d' % mezz]['overload_dmfd_%s' % (['a', 'b'][mod-1])]

                modhk.module_number = mod

                for chan in range(1, dat['config']['dmfd_channels_per_wire'] + 1):
                    chanhk = HkChannelInfo()
                    chanhk.channel_number = chan
                    chanhk.dan_gain = dat['dan']['gain'][wire-1][chan-1]
                    chanhk.dan_accumulator_enable = dat['dan']['accumulator_enable'][wire-1][chan-1]
                    chanhk.dan_feedback_enable = dat['dan']['feedback_enable'][wire-1][chan-1]
                    chanhk.dan_streaming_enable = dat['dan']['streaming_enable'][wire-1][chan-1]
                    chanhk.dan_railed = dat['dan']['accumulator_railed'][wire-1][chan-1]

                    # Channel params: (freq, phase, amp)
                    carrier = dat['dmfs']['carrier'][wire-1][chan-1]
                    nuller = dat['dmfs']['nuller'][wire-1][chan-1]
                    # Demod params: (freq, phase)
                    demod = dat['dmfd']['demod'][wire-1][chan-1]
                    chanhk.carrier_amplitude = carrier[2]/(2.**23 - 1.)
                    chanhk.nuller_amplitude = nuller[2]/(2.**23 - 1.)
                    chanhk.carrier_frequency = carrier[0]*25.e6/(2.**32-1)*core.G3Units.Hz
                    chanhk.demod_frequency = demod[0]*25.e6/(2.**32-1)*core.G3Units.Hz

                    modhk.channels[chan] = chanhk
                mezzhk.modules[mod] = modhk
            boardhk.mezz[mezz] = mezzhk
    
        return boardhk

    @classmethod
    def HousekeepingForBoard(cls, hostname):
        '''
        Return HkBoardInfo object from a board. Useful for debugging.
        '''
        t = TuberClient(hostname)
        jsonblob = t.CallMethod('Dfmux', ['classic_dfmux_dump', 'classic_housekeeping'])
        data = dict(jsonblob[0]['result'])
        data.update(jsonblob[1]['result'])

        return cls.HousekeepingFromJSON(data)

