from spt3g import core
from spt3g.dfmux import DfMuxHousekeepingMap, HkBoardInfo, HkMezzanineInfo, HkModuleInfo, HkChannelInfo

from spt3g.dfmux.IceboardConversions import convert_TF
from .TuberClient import TuberClient
import socket, struct, time
import numpy

@core.indexmod
class HousekeepingConsumer(object):
    '''
    Collect housekeeping data from the mux boards defined in the wiring map.
    Will add a key called 'DfMuxHousekeeping' to any housekeeping frame that
    goes by containing the data as of the arrival of the housekeeping frame.
    Use in conjunction with a dfmux.PeriodicHousekeepingCollector to get
    data at fixed intervals.
   
    If collecting real-time data, you may want to set subprocess=True when
    adding this module.
    '''
    def __init__(self):
        self.tuber = None
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Wiring:
            wmap = frame['WiringMap']
            self.board_ips = set([m.board_ip for m in wmap.values()])
            self.tuber = [(ip, TuberClient(socket.inet_ntoa(struct.pack('i', ip)), timeout=5.0)) for ip in self.board_ips]

        if frame.type == core.G3FrameType.Housekeeping:
            if self.tuber is None:
                self.tuber = [(ip, TuberClient(socket.inet_ntoa(struct.pack('i', ip)), timeout=5.0)) for ip in self.board_ips]

            hkdata = DfMuxHousekeepingMap()
            try:
                for board,board_tuber in self.tuber:
                    board_tuber.CallMethod('Dfmux', '_dump_housekeeping', False)
                    time.sleep(0.02) # Stagger return data transfer a little to
                                     # avoid overloading the network on the return

                for board,board_tuber in self.tuber:
                    dat = board_tuber.GetReply()[0]['result']

                    boardhk = self.HousekeepingFromJSON(dat)
                    hkdata[int(boardhk.serial)] = boardhk
                frame['DfMuxHousekeeping'] = hkdata
            except socket.timeout:
                core.log_error('Timeout collecting housekeeping data from mux boards. Dropping housekeeping sample', unit='HousekeepingConsumer')
                return []
            except Exception as e:
                core.log_error('Error (%s) collecting housekeeping data from mux boards. Dropping housekeeping sample' % e, unit='HousekeepingConsumer')
                return []

    @classmethod
    def HousekeepingFromJSON(cls, dat):
        '''
        Build HKBoardInfo object from a JSON blob returned by the
        _dump_housekeeping call
        '''
        # Board-global quantities
        boardhk = HkBoardInfo()

        year = dat['timestamp']['y']
        if year == 0:
            # It probably isn't 1900
            systime = time.gmtime()
            # Check for New Year's, assuming no more than 24 hours clock slew
            if dat['timestamp']['d'] == 1 and systime.tm_yday >= 365:
                year = systime.tm_year + 1
            elif dat['timestamp']['d'] >= 365 and systime.tm_yday == 1:
                year = systime.tm_year - 1
            else:
                year = systime.tm_year
            year -= 2000
        boardhk.timestamp = core.G3Time(y=year,d=dat['timestamp']['d'],h=dat['timestamp']['h'],m=dat['timestamp']['m'],s=dat['timestamp']['s'],ss=dat['timestamp']['ss'])

        boardhk.timestamp_port = str(dat['timestamp_port'])
        boardhk.serial = str(dat['serial'])
        boardhk.fir_stage = dat['fir_stage']
        for i in dat['currents'].items():
            boardhk.currents[str(i[0])] = i[1]
        for i in dat['voltages'].items():
            boardhk.voltages[str(i[0])] = i[1]
        for i in dat['temperatures'].items():
            boardhk.temperatures[str(i[0])] = i[1]

        # Mezzanines
        for n, mezz in enumerate(dat['mezzanines']):
            mezzhk = HkMezzanineInfo()
            mezzhk.present = mezz['present']
            mezzhk.power = mezz['power']
            if mezzhk.present:
                mezzhk.serial = str(mezz['ipmi']['product']['serial_number'])
                mezzhk.part_number = str(mezz['ipmi']['product']['part_number'])
                mezzhk.revision = str(mezz['ipmi']['product']['version_number'])
                for i in mezz['currents'].items():
                    mezzhk.currents[str(i[0])] = i[1]
                for i in mezz['voltages'].items():
                    mezzhk.voltages[str(i[0])] = i[1]

            if mezzhk.present and mezzhk.power:
                mezzhk.temperature = mezz['temperature']
                mezzhk.squid_heater = 0.0 #mezz['squid_heater']  # AJA/KTS: temporary kludge, not in housekeeping tuber
                mezzhk.squid_controller_power = False #mezz['squid_controller_power']  # AJA/KTS: temporary kludge, not in housekeeping tuber
                mezzhk.squid_controller_temperature = 0.0 #mezz['squid_controller_temperature']  # AJA/KTS: temporary kludge, not in housekeeping tuber

            # Modules
            for m, mod in enumerate(mezz['modules']):
                modhk = HkModuleInfo()
                modhk.routing = str(mod['routing'][0])
                modhk.module_number = m+1
                if mezzhk.present and mezzhk.power:
                    if 'gains' in mod:
                        modhk.carrier_gain = mod['gains']['carrier']
                        modhk.nuller_gain = mod['gains']['nuller']
                    if 'overload' in mod:
                        modhk.carrier_railed = mod['overload']['carrier']
                        modhk.nuller_railed = mod['overload']['nuller']
                        modhk.demod_railed = mod['overload']['demod']
                    if 'squid_current_bias' in mod:
                        modhk.squid_current_bias = mod['squid_current_bias']
                    if 'squid_flux_bias' in mod:
                        modhk.squid_current_bias = mod['squid_flux_bias']
                    if 'squid_feedback' in mod:
                        modhk.squid_feedback = str(mod['squid_feedback'])

                if 'squid_tuning' in mod and mod['squid_tuning'] is not None:
                    modhk.squid_state = str(mod['squid_tuning']['state'])
                    modhk.squid_transimpedance = mod['squid_tuning']['transimpedance'] if mod['squid_tuning']['transimpedance'] is not None else numpy.nan
                    modhk.squid_p2p = mod['squid_tuning']['p2p'] if mod['squid_tuning']['p2p'] is not None else numpy.nan

                for k, chan in enumerate(mod['channels']):
                    chanhk = HkChannelInfo()
                    chanhk.channel_number = k+1
                    chanhk.carrier_amplitude = chan['carrier_amplitude']
                    chanhk.nuller_amplitude = chan['nuller_amplitude']
                    chanhk.carrier_frequency = chan['carrier_frequency']*core.G3Units.Hz
                    chanhk.demod_frequency = chan['demod_frequency']*core.G3Units.Hz
                    chanhk.dan_gain = chan['dan_gain']
                    chanhk.dan_accumulator_enable = chan['dan_accumulator_enable']
                    chanhk.dan_feedback_enable = chan['dan_feedback_enable']
                    chanhk.dan_streaming_enable = chan['dan_streaming_enable']
                    if 'dan_railed' in chan:
                        chanhk.dan_railed = chan['dan_railed']
                    if 'tuning' in chan and chan['tuning'] is not None:
                        chanhk.state = str(chan['tuning']['state'])
                        if ('rlatched' in chan['tuning'] and 
                            chan['tuning']['rlatched'] != None):
                            chanhk.rlatched = chan['tuning']['rlatched']
                        if ('rnormal' in chan['tuning'] and 
                            chan['tuning']['rnormal'] != None):
                            chanhk.rnormal = chan['tuning']['rnormal']
                        if ('rfrac_achieved' in chan['tuning'] and 
                            chan['tuning']['rfrac_achieved'] != None):
                            chanhk.rfrac_achieved = chan['tuning']['rfrac_achieved']

                    #calculate the resistance conversion factor
                    if 'gains' in mod:
                        V = convert_TF(modhk.carrier_gain, 
                                       target='carrier', unit='NORMALIZED', 
                                       frequency = chanhk.carrier_frequency / core.G3Units.Hz)
                        I = convert_TF(modhk.nuller_gain, 
                                       target='nuller', unit='RAW', 
                                       frequency = chanhk.carrier_frequency / core.G3Units.Hz)
                        chanhk.res_conversion_factor = chanhk.carrier_amplitude * V/I
                        #  R = (V * Vconv) / (I  * Iconv)

                    modhk.channels[k+1] = chanhk
                mezzhk.modules[m+1] = modhk
            boardhk.mezz[n+1] = mezzhk
    
        return boardhk

    @classmethod
    def HousekeepingForBoard(cls, hostname):
        '''
        Return HkBoardInfo object from a board. Useful for debugging.
        '''
        t = TuberClient(hostname)
        data = t.CallMethod('Dfmux', '_dump_housekeeping')
        return cls.HousekeepingFromJSON(data[0]['result'])

@core.indexmod
class PeriodicHousekeepingCollector(object):
    '''Inserts housekeeping frames every N timepoints.'''
    def __init__(self, N=15200):
        self.N = N
        self.count = 0
    def __call__(self, frame):
        ret = []
        if frame.type == core.G3FrameType.Timepoint:
            if self.count % self.N == 0:
                ret.append(core.G3Frame(core.G3FrameType.Housekeeping))
            self.count += 1
            ret.append(frame)
            return ret

@core.usefulfunc
def HousekeepingForBolo(hkmap, wiringmap, bolo, all_hk=False):
    '''
    Obtain the channel housekeeping information for a bolometer named "bolo"
    using the passed housekeeping and wiring maps.

    If all_hk is True, returns a tuple of the (board, mezz, module, channel)
    HK data instead of just the channel.
    '''

    wiringentry = wiringmap[bolo]
    board = hkmap[wiringentry.board_serial]
    mezz = board.mezz[(wiringentry.module // len(board.mezz[1].modules)) + 1]
    mod = mezz.modules[(wiringentry.module % len(mezz.modules)) + 1]
    chan = mod.channels[wiringentry.channel + 1]
    if all_hk:
        return (board, mezz, mod, chan)
    else:
        return chan
