from spt3g import core
from spt3g.dfmux import DfMuxHousekeepingMap, HkBoardInfo, HkMezzanineInfo, HkModuleInfo, HkChannelInfo, DfMuxWiringMap, DfMuxChannelMapping

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

    Also emits a Wiring frame from the housekeeping data.  This requires a
    recent (as of November 2018) version of pydfmux in order to read mapped
    channel names from each board.  If ignore_wiring=False, assumes that
    the wiring map is handled by a separate process.

    If collecting real-time data, you may want to set subprocess=True when
    adding this module.
    '''
    def __init__(self, ignore_wiring=False):
        self.tuber = {}
        self.board_map = {}
        self.board_serials = []
        self.buffer = []
        self.hwmf = None
        self.ignore_wiring = ignore_wiring

    def map_boards(self):
        '''
        Cache IP address, crate serial and slot number and
        create a TuberClient for each IceBoard.
        '''
        boards = [b for b in self.board_serials if b not in self.board_map]
        if not boards:
            return

        ips = {}
        crates = {}
        slots = {}

        # create tubers for each missing board
        for board in boards:
            ip = socket.inet_aton(socket.gethostbyname('iceboard{}.local'.format(board)))
            ips[board] = struct.unpack("i", ip)[0]
            self.tuber[board] = TuberClient(socket.inet_ntoa(ip), timeout=5.0)

        # get crate ID for each missing board
        for board in boards:
            self.tuber[board].CallMethod('Dfmux', '_get_backplane_serial', False)
        for board in boards:
            reply = self.tuber[board].GetReply()[0]
            if reply.get('error', None):
                crates[board] = -1
            else:
                crates[board] = int(reply['result'])

        # get slot ID for each missing board
        for board in boards:
            if crates[board] >= 0:
                self.tuber[board].CallMethod('Dfmux', 'get_backplane_slot', False)
        for board in boards:
            if crates[board] < 0:
                slots[board] = -1
            else:
                reply = self.tuber[board].GetReply()[0]
                if reply.get('error', None):
                    slots[board] = -1
                else:
                    slots[board] = int(reply['result'])

        # store
        for board in boards:
            self.board_map[board] = (ips[board], crates[board], slots[board])

    def __call__(self, frame):

        # If boards have been mapped already, then process every frame
        # as received.
        if len(self.board_map):
            return self.ProcessBuffered(frame)

        # Otherwise, process at least one Timepoint frame first to gather
        # a list of IceBoards for which housekeeping data are required,
        # Then process any buffered frames after Timepoint frame,
        # but emit before the Timepoint
        if frame.type == core.G3FrameType.Timepoint:
            tp_frames = self.ProcessBuffered(frame)
            frames = sum([self.ProcessBuffered(fr) for fr in self.buffer], [])
            frames += tp_frames
            self.buffer = []
            return frames

        # Buffer until a Timepoint frame is received
        self.buffer.append(frame)
        return False

    def ProcessBuffered(self, frame):
        '''
        Process frames in the buffered order.  Returned value should
        always be a list of frames, possibly empty.
        '''
        if frame.type == core.G3FrameType.Timepoint:
            self.board_serials = [('%04d' % k) for k in frame['DfMux'].keys()]
            return [frame]

        if frame.type == core.G3FrameType.Wiring:
            if self.ignore_wiring:
                core.log_warn(
                    "Received a wiring frame, which may be inconsistent with "
                    "board housekeeping data.  Store mapped channel names on the "
                    "board using updated pydfmux/hidfmux software.",
                    unit="HousekeepingConsumer",
                )
            else:
                core.log_fatal(
                    "Received spurious wiring frame.  Do not use "
                    "PyDfMuxHardwareMapInjector with the HousekeepingConsumer.  "
                    "You may update pydfmux/hidfmux to a newer version that "
                    "stores mapped channel names on the boards, and rearrange "
                    "your data acquisition script.",
                    unit="HousekeepingConsumer",
                )

        if frame.type == core.G3FrameType.Housekeeping:
            self.map_boards()

            hwm = DfMuxWiringMap()
            hkdata = DfMuxHousekeepingMap()
            try:
                for board in self.board_serials:
                    self.tuber[board].CallMethod('Dfmux', '_dump_housekeeping', False)
                    time.sleep(0.02) # Stagger return data transfer a little to
                                     # avoid overloading the network on the return

                found = False
                ismkid = False
                for board in self.board_serials:
                    dat = self.tuber[board].GetReply()[0]['result']

                    boardhk = self.HousekeepingFromJSON(dat)
                    if boardhk.firmware_name:
                        ismkid = ismkid or "mkid" in boardhk.firmware_name.lower()
                    hkdata[int(boardhk.serial)] = boardhk

                    if self.ignore_wiring:
                        continue

                    ip, crate, slot = self.board_map[board]
                    boardw = self.WiringFromJSON(dat, ip, crate, slot)
                    for key in boardw.keys():
                        hwm[key] = boardw[key]
                        found = True

                if not found and not self.ignore_wiring:
                    core.log_fatal("No mapped channels found on any IceBoards. "
                                   "You may need to update pydfmux to a newer version of pydfmux "
                                   "that stores mapped channel names on the boards, and reload "
                                   "the hardware map.",
                                   unit='HousekeepingConsumer')

                frame['DfMuxHousekeeping'] = hkdata

                if self.ignore_wiring:
                    return [frame]

                hwmf = core.G3Frame(core.G3FrameType.Wiring)
                hwmf['WiringMap'] = hwm
                hwmf['ReadoutSystem'] = 'RF-ICE' if ismkid else 'ICE'

                if self.hwmf is None:
                    self.hwmf = hwmf
                    # If this is the first time the consumer is triggered, make sure
                    # a Housekeeping and WiringMap frame are issued before any
                    # Timepoint frames
                    frames = [hwmf, frame]
                    return frames
                else:
                    # Compare wiring maps and re-issue frame if anything has changed
                    old_hwm = self.hwmf['WiringMap']
                    if (set(hwm.keys()) ^ set(old_hwm.keys())):
                        self.hwmf = hwmf
                        return [hwmf, frame]
                    for k in hwm.keys():
                        try:
                            if vars(hwm[str(k)]) != vars(old_hwm[str(k)]):
                                self.hwmf = hwmf
                                return [hwmf, frame]
                        except:
                            core.log_error("Invalid HWM key %r" % k)

                # If we get here then the wiring map hasn't changed,
                # so return the populated Housekeeping frame as it is
                return [frame]

            except socket.timeout:
                core.log_error('Timeout collecting housekeeping data from mux boards. Dropping housekeeping sample', unit='HousekeepingConsumer')
                return []
            except Exception as e:
                core.log_error('Error (%s) collecting housekeeping data from mux boards. Dropping housekeeping sample' % e, unit='HousekeepingConsumer')
                return []

        return [frame]

    @classmethod
    def HousekeepingFromJSON(cls, dat):
        '''
        Build HKBoardInfo object from a JSON blob returned by the
        _dump_housekeeping call
        '''
        # Board-global quantities
        boardhk = HkBoardInfo()
        if 'is128x' in dat:
            boardhk.is128x = dat['is128x']
        else:
            boardhk.is128x = False
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
        if 'firmware_name' in dat:
            boardhk.firmware_name = str(dat['firmware_name'])
            boardhk.firmware_version = str(dat['firmware_version'])
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
                if 'ipmi' in mezz:
                    mezzhk.serial = str(mezz['ipmi']['product']['serial_number'])
                    mezzhk.part_number = str(mezz['ipmi']['product']['part_number'])
                    mezzhk.revision = str(mezz['ipmi']['product']['version_number'])
                if 'currents' in mezz:
                    for i in mezz['currents'].items():
                        mezzhk.currents[str(i[0])] = i[1]
                if 'voltages' in mezz:
                    for i in mezz['voltages'].items():
                        mezzhk.voltages[str(i[0])] = i[1]

            if mezzhk.present and mezzhk.power:
                if 'temperature' in mezz:
                    mezzhk.temperature = mezz['temperature']
                # these parameters are not in the 64x housekeeping tuber
                mezzhk.squid_heater = mezz.get('squid_heater', 0.0)
                mezzhk.squid_controller_power = mezz.get('squid_controller_power', False)
                mezzhk.squid_controller_temperature = mezz.get('squid_controller_temperature', 0.0)

            # Modules
            for m, mod in enumerate(mezz['modules']):
                modhk = HkModuleInfo()
                modhk.routing_type = str(mod['routing'][0])
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
                    if 'nco_frequency' in mod:
                        modhk.nco_frequency = mod['nco_frequency'] * core.G3Units.Hz

                if 'squid_tuning' in mod and mod['squid_tuning'] is not None:
                    modhk.squid_state = str(mod['squid_tuning']['state'])
                    modhk.squid_transimpedance = mod['squid_tuning']['transimpedance'] if mod['squid_tuning']['transimpedance'] is not None else numpy.nan
                    modhk.squid_p2p = mod['squid_tuning']['p2p'] if mod['squid_tuning']['p2p'] is not None else numpy.nan

                for k, chan in enumerate(mod['channels']):
                    chanhk = HkChannelInfo()
                    chanhk.channel_number = k+1
                    chanhk.carrier_amplitude = chan['carrier_amplitude']
                    chanhk.nuller_amplitude = chan['nuller_amplitude']
                    if 'dan_gain' in chan:
                        chanhk.dan_gain = chan['dan_gain']
                    if 'dan_streaming_enable' in chan:
                        chanhk.dan_streaming_enable = chan['dan_streaming_enable']
                    if 'frequency' in chan:
                        chanhk.carrier_frequency = chan['frequency'] * core.G3Units.Hz
                        chanhk.demod_frequency = chan['frequency'] * core.G3Units.Hz
                    else:
                        chanhk.carrier_frequency = chan['carrier_frequency'] * core.G3Units.Hz
                        chanhk.demod_frequency = chan['demod_frequency'] * core.G3Units.Hz
                    if 'carrier_phase' in chan:
                        chanhk.carrier_phase = chan['carrier_phase'] * core.G3Units.deg
                        chanhk.nuller_phase = chan['nuller_phase'] * core.G3Units.deg
                        chanhk.demod_phase = chan['demod_phase'] * core.G3Units.deg
                    if 'dan_accumulator_enable' in chan:
                        chanhk.dan_accumulator_enable = chan['dan_accumulator_enable']
                    if 'dan_feedback_enable' in chan:
                        chanhk.dan_feedback_enable = chan['dan_feedback_enable']
                    if 'dan_railed' in chan:
                        chanhk.dan_railed = chan['dan_railed']
                    if 'tuning' in chan and chan['tuning'] is not None:
                        chanhk.state = str(chan['tuning']['state'])
                        if ('rlatched' in chan['tuning'] and 
                            chan['tuning']['rlatched'] is not None):
                            chanhk.rlatched = chan['tuning']['rlatched']
                        if ('rnormal' in chan['tuning'] and 
                            chan['tuning']['rnormal'] is not None):
                            chanhk.rnormal = chan['tuning']['rnormal']
                        if ('rfrac_achieved' in chan['tuning'] and 
                            chan['tuning']['rfrac_achieved'] is not None):
                            chanhk.rfrac_achieved = chan['tuning']['rfrac_achieved']
                        if ('loopgain' in chan['tuning'] and
                            chan['tuning']['loopgain'] is not None):
                            chanhk.loopgain = chan['tuning']['loopgain']

                    modhk.channels[k+1] = chanhk
                mezzhk.modules[m+1] = modhk
            boardhk.mezz[n+1] = mezzhk
    
        return boardhk

    @classmethod
    def WiringFromJSON(cls, dat, ip, crate, slot):
        '''
        Build WiringMap object from a JSON blob returned by the
        _dump_housekeeping call
        '''
        found = False
        hwm = DfMuxWiringMap()
        serial = int(dat['serial'])
        for imezz, mezz in enumerate(dat['mezzanines']):
            for imod, mod in enumerate(mezz['modules']):
                module = imod + len(mezz['modules']) * imezz
                for ichan, chan in enumerate(mod['channels']):
                    name = (chan.get('tuning', {}) or {}).get('name', None)
                    if not name:
                        continue
                    mapping = DfMuxChannelMapping()
                    mapping.board_ip = ip
                    mapping.board_serial = serial
                    mapping.board_slot = slot
                    mapping.crate_serial = crate
                    mapping.module = module
                    mapping.channel = ichan
                    try:
                        name = str(name)
                        hwm[name] = mapping
                        found = True
                    except:
                        core.log_error("Invalid channel name %r" % (name))
        if not found:
            core.log_error("No mapped channels found on iceboard%04d. "
                           "You may need to update pydfmux to a newer version of pydfmux "
                           "that stores mapped channel names on the boards, and reload "
                           "the hardware map." % (serial),
                           unit='HousekeepingConsumer')
        return hwm

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
