import struct, socket, errno, numpy, time, threading
from spt3g import core, dfmux


class PagerWatchdog(object):
    '''
    Module that sends a watchdog (ping) message to the GCP pager when the parent
    process is running successfully.  Modify the `data_valid` method for
    particular use cases, and call the `run` method periodically in your
    application.
    '''

    host = 'sptnet.spt'
    port = 50040
    timeout = 20

    def __init__(self, name, interval=600, sim=False):
        self.name = name.lower()
        self.unit = '{}Watchdog'.format(name.capitalize())
        self.interval = interval
        self.sim = sim
        self.last_ping = None

        # ping on startup
        self.thread = threading.Thread(target=self.ping)
        self.thread.start()

    def ping(self):
        """
        Send a watchdog ping message to the GCP pager process.  This method is
        called by the `run` method at regular intervals whenever the
        `data_valid` method returns True.
        """
        try:
            if not self.sim:
                sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                sock.settimeout(self.timeout)
                sock.connect((self.host, self.port))
                sock.send('watchdog {}'.format(self.name).encode())
                resp = sock.recv(4096)
                if resp:
                    core.log_debug(
                        'Sent watchdog ping, got response {}'.format(resp.decode()),
                        unit=self.unit,
                    )
                sock.close()
        except Exception as e:
            core.log_error('Error sending watchdog ping: {}'.format(e), unit=self.unit)
            # try again in ten seconds
            self.last_ping = time.time() - self.interval + 10
        else:
            core.log_info('Sent watchdog ping', unit=self.unit)
            self.last_ping = time.time()

    def data_valid(self, *args, **kwargs):
        """
        Returns True if the watchdog should ping, otherwise False.  Arguments
        are passed to this method from the `run` method.
        """
        raise NotImplementedError

    def run(self, *args, **kwargs):
        """
        When called, issues a watchdog ping message if the interval time has passed, and
        the `data_valid` method returns True.  All input arguments are passed to the
        `data_valid` method for validation.
        """

        # only ping if ready
        if not self.data_valid(*args, **kwargs):
            return

        # only ping if another ping isn't already running
        if self.thread is not None:
            if not self.thread.is_alive():
                del self.thread
                self.thread = None
            else:
                return

        # only ping on the appropriate interval
        now = time.time()
        if self.last_ping and (now - self.last_ping < self.interval):
            return

        # ping
        self.thread = threading.Thread(target=self.ping)
        self.thread.start()


@core.indexmod
class DAQWatchdog(PagerWatchdog):
    """
    Watchdog that issues a ping to the GCP pager when the DAQ is operational.
    """

    def __init__(self, calibrator=False, interval=600, sim=False):
        """
        Arguments
        ---------
        calibrator : bool
            If True, ensure that the calibrator is also running successfully
            before sending a ping.
        """
        super(DAQWatchdog, self).__init__('DAQ', interval=interval, sim=sim)

        self.last_missing = None
        self.boards_missing = 0
        self.last_log_boards = None
        self.calibrator = calibrator
        self.last_log_cal = None

    def data_valid(self, frame):
        """
        Check the incoming frame for completeness.

         * Ensure that all modules in the listed iceboards are reporting.
         * If `calibrator` is True, ensure that the calibrator sync signal is in the frame.
        """

        # always ready in sim mode
        if self.sim:
            return True

        # only ping on Timepoint frames
        if 'DfMux' not in frame:
            return False

        now = time.time()
        retval = True

        # only ping if all expected modules are present
        data = frame['DfMux']
        boards_expected = len(data.keys())
        boards_complete = sum([v.nmodules > 0 and v.Complete() for v in data.values()])
        boards_missing = boards_expected - boards_complete
        if boards_missing:
            if not self.last_log_boards or boards_missing != self.boards_missing:
                # log loss or change in missing count
                core.log_error(
                    "Missing data from {} boards in DAQ data stream".format(boards_missing),
                    unit=self.unit,
                )
                self.last_log_boards = now
                self.boards_missing = boards_missing
            self.last_missing = now
            retval = False
        elif self.last_log_boards and now - self.last_missing > 10:
            # log recovery
            core.log_notice("All boards recovered in DAQ data stream", unit=self.unit)
            self.last_log_boards = None

        # only ping if the calibrator sync signal is present
        if self.calibrator and 'CalibratorOn' not in frame:
            if not self.last_log_cal:
                # log loss
                core.log_error(
                    "Missing calibrator signal in DAQ data stream",
                    unit=self.unit,
                )
                self.last_log_cal = now
            self.last_missing = now
            retval = False
        elif self.calibrator and self.last_log_cal and now - self.last_missing > 10:
            # log recovery
            core.log_notice(
                "Calibrator signal recovered in DAQ data stream", unit=self.unit
            )
            self.last_log_cal = None

        # only ping if normal data acquisition has been going for a bit
        if self.last_missing and now - self.last_missing < 10:
            retval = False

        return retval

    def __call__(self, frame):
        self.run(frame)


@core.indexmod
class GCPSignalledHousekeeping(object):
    '''
    Module that collects housekeeping data when connected to. If
    collect_on_start is True (the default), injects an HK frame
    unconditionally at startup.
    '''
    def __init__(self, port=50011, collect_on_start=True):
        self.socket = socket.socket()
        self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)

        self.socket.bind(('', port))
        self.socket.listen(5)
        self.socket.setblocking(False)

        self.collect_on_start = collect_on_start
        self.first_frame = True

    def __del__(self):
        self.socket.close()

    def __call__(self, frame):
        # Only try on timepoint frames
        if frame.type != core.G3FrameType.Timepoint:
            return

        if self.collect_on_start and self.first_frame:
            self.first_frame = False
            return [core.G3Frame(core.G3FrameType.Housekeeping), frame]

        # Check for new connections
        try:
            s, origin_ip = self.socket.accept()
        except socket.error as e:
            if e.errno != errno.EAGAIN and e.errno != errno.EWOULDBLOCK:
                raise
            return
        core.log_debug('Accepted housekeeping collection signal from %s:%d' % origin_ip,
                       unit='GCPSignalledHousekeeping')
        s.close()

        return [core.G3Frame(core.G3FrameType.Housekeeping), frame]

@core.indexmod
class GCPHousekeepingTee(object):
    '''
    Module that serves housekeeping information to GCP when asked. If a key
    named "DataOK" exists in the housekeeping frames, will also transmit
    that information to GCP for paging purposes.
    '''
    def __init__(self, port=50010, verbose=False):

        # make some noise at startup
        core.log_info("Initialize gcp.GCPHousekeepingTee on port %d" % port,
                      unit='GCPHousekeepingTee')

        self.hkblob = self.PackHKToGCP(dfmux.DfMuxHousekeepingMap())
        self.socket = socket.socket()
        self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)

        self.socket.bind(('', port))
        self.socket.listen(25)
        self.socket.setblocking(False)

        self.verbose = verbose  # flag for printing debugging statements

    def __del__(self):
        # Clear any pending connections. No one is getting anything now.
        # This works around some misfeatures in the Linux kernel that
        # do not occur in other, better socket implementations.
        while True:
            try:
                s, origin_ip = self.socket.accept()
            except socket.error as e:
                if e.errno != errno.EAGAIN and e.errno != errno.EWOULDBLOCK:
                    raise
                break

            s.close()
        self.socket.close()

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Housekeeping:
            dataok = True
            if 'DataOK' in frame:
                dataok = frame['DataOK']
            self.hkblob = self.PackHKToGCP(
                frame['DfMuxHousekeeping'], dataok=dataok,
                verbose=self.verbose)

            # Check for new connections, send any interested
            # parties the same data
            cxs = []
            while True:
                try:
                    s, origin_ip = self.socket.accept()
                except socket.error as e:
                    if e.errno != errno.EAGAIN and \
                       e.errno != errno.EWOULDBLOCK:
                        raise
                    break

                core.log_debug('Accepted connection from %s:%d' % origin_ip,
                               unit='GCPHousekeepingTee')
                cxs.append(s)

            for s in cxs:
                s.setblocking(True)
                s.sendall(self.hkblob)
                s.close()


    @staticmethod
    def PackHKToGCP(hk, dataok=True, verbose=False):
        if verbose:
            core.log_debug('gcp.GCPHousekeepingTee.PackHKToGCP(hk)', unit='GCPHousekeepingTee')
        buf = struct.pack('<?I', dataok, len(hk))

        # See HkDataStruct in GCP
        for ip, board in hk.items():

            # if verbose mode, print a few registers for debugging
            if verbose:
                core.log_debug("%d, %10.6f, %10.6f, %10.6f, %10.6f, %10.6f" %
                               (ip, board.temperatures['MOTHERBOARD_TEMPERATURE_ARM'],
                                board.temperatures['MOTHERBOARD_TEMPERATURE_FPGA'],
                                board.temperatures['MOTHERBOARD_TEMPERATURE_FPGA_DIE'],
                                board.temperatures['MOTHERBOARD_TEMPERATURE_PHY'],
                                board.temperatures['MOTHERBOARD_TEMPERATURE_POWER']),
                               unit='GCPHousekeepingTee')

            buf += struct.pack('<fffff',
                               board.temperatures['MOTHERBOARD_TEMPERATURE_ARM'],
                               board.temperatures['MOTHERBOARD_TEMPERATURE_FPGA'],
                               board.temperatures['MOTHERBOARD_TEMPERATURE_FPGA_DIE'],
                               board.temperatures['MOTHERBOARD_TEMPERATURE_PHY'],
                               board.temperatures['MOTHERBOARD_TEMPERATURE_POWER'])
            buf += struct.pack('<fffffffff',
                               board.voltages['MOTHERBOARD_RAIL_VADJ'],
                               board.voltages['MOTHERBOARD_RAIL_VCC12V0'],
                               board.voltages['MOTHERBOARD_RAIL_VCC1V0'],
                               board.voltages['MOTHERBOARD_RAIL_VCC1V0_GTX'],
                               board.voltages['MOTHERBOARD_RAIL_VCC1V2'],
                               board.voltages['MOTHERBOARD_RAIL_VCC1V5'],
                               board.voltages['MOTHERBOARD_RAIL_VCC1V8'],
                               board.voltages['MOTHERBOARD_RAIL_VCC3V3'],
                               board.voltages['MOTHERBOARD_RAIL_VCC5V5'])
            buf += struct.pack('<fffffffff',
                               board.currents['MOTHERBOARD_RAIL_VADJ'],
                               board.currents['MOTHERBOARD_RAIL_VCC12V0'],
                               board.currents['MOTHERBOARD_RAIL_VCC1V0'],
                               board.currents['MOTHERBOARD_RAIL_VCC1V0_GTX'],
                               board.currents['MOTHERBOARD_RAIL_VCC1V2'],
                               board.currents['MOTHERBOARD_RAIL_VCC1V5'],
                               board.currents['MOTHERBOARD_RAIL_VCC1V8'],
                               board.currents['MOTHERBOARD_RAIL_VCC3V3'],
                               board.currents['MOTHERBOARD_RAIL_VCC5V5'])
            buf += struct.pack('255s', ('iceboard' + board.serial).encode())
            for i in [1,2]:
                buf += struct.pack('<?',
                                   board.mezz[i].present)
                buf += struct.pack('<?',
                                   board.mezz[i].power)
                buf += struct.pack('<?',
                                   board.mezz[i].squid_controller_power)
                buf += struct.pack('<fff',
                                   board.mezz[i].voltages['MEZZANINE_RAIL_VADJ'],
                                   board.mezz[i].voltages['MEZZANINE_RAIL_VCC12V0'],
                                   board.mezz[i].voltages['MEZZANINE_RAIL_VCC3V3'])
                buf += struct.pack('<fff',
                                   board.mezz[i].currents['MEZZANINE_RAIL_VADJ'],
                                   board.mezz[i].currents['MEZZANINE_RAIL_VCC12V0'],
                                   board.mezz[i].currents['MEZZANINE_RAIL_VCC3V3'])
                buf += struct.pack('<fff',
                                   board.mezz[i].temperature,
                                   board.mezz[i].squid_controller_temperature,
                                   board.mezz[i].squid_heater)


        # Prefix with total message length
        buf = struct.pack('!q', len(buf)) + buf

        return buf

@core.indexmod
class GCPBoloDataTee(object):
    '''
    Module that serves bolometer data to GCP when asked. Once a second,
    will serve the data from the previous second of bolometer data.

    If a boolean key appears in the timepoint frames named "DataOK",
    this will be sent to GCP as a data quality indicator for paging.
    '''
    def __init__(self, port=50020, verbose=False, bolometers=[]):
        '''
        Send data from the given list of bolometer logical IDs to the GCP.
        '''

        core.log_info('Listening for requests from GCP on port %d' % port, unit='GCPBoloDataTee')
        core.log_info('Selected bolometers: %s' % bolometers, unit='GCPBoloDataTee')

        self.bololist = bolometers

        self.socket = socket.socket()
        self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)

        self.socket.bind(('', port))
        self.socket.listen(5)
        self.socket.setblocking(False)
        self.data = {}

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Wiring:
            self.wmap = frame['WiringMap']

        if frame.type == core.G3FrameType.Timepoint:
            sec = int(frame['EventHeader'].time/core.G3Units.s)

            # Add data from this sample to the cache for this calendar second, replacing any missing detectors with -1
            if sec not in self.data:
                self.data[sec] = {b: [] for b in self.bololist}
                self.data[sec]['DataOK'] = []
            d = self.data[sec]
            for b in self.bololist:
                w = self.wmap[b]
                board = frame['DfMux'][w.board_serial]
                if board.nblocks > 1:
                    mod_idx = w.module * board.nblocks + w.channel // board.nchannels
                    chan_idx = w.channel % board.nchannels
                else:
                    mod_idx = w.module
                    chan_idx = w.channel
                try:
                    d[b].append(board[mod_idx][chan_idx])
                except KeyError:
                    d[b].append(-1)
            if 'DataOK' in frame:
                self.data[sec]['DataOK'].append(bool(frame['DataOK']))
            else:
                self.data[sec]['DataOK'].append(True)

            # Toss ancient data: we keep the last second (complete)
            # for GCP, plus the second we are currently accumulating
            if len(self.data) > 2:
                keys = list(self.data.keys())
                keys.sort()
                for k in keys[:-2]:
                    del self.data[k]

        # Check for new connections once we have a buffer
        if len(self.data) == 2:
            try:
                s, origin_ip = self.socket.accept()
            except socket.error as e:
                if e.errno != errno.EAGAIN and e.errno != errno.EWOULDBLOCK:
                    raise
                return

            core.log_debug('Accepted connection from %s:%d' % origin_ip, unit='GCPBoloDataTee')
            s.setblocking(True)
            keys = list(self.data.keys())
            keys.sort()
            s.sendall(self.PackForGCP(self.data[keys[0]]))
            s.close()

            # Delete data once enqueued
            del self.data[keys[0]]


    @staticmethod
    def PackForGCP(data):
        # Input data: dict of bolo names to samples (can be a
        # G3TimestreamMap, in principle)
        #
        # On-wire format:
        # U64 Length of buffer (big-endian)
        # U8 Data Quality Indicator: True (1) = good
        # U32 Number of detectors in list
        # U32 Number of samples in the last second
        # N copies of:
        # - 16 byte character string with detector name
        # - N_sample 32-bit signed integers with data
        #
        # All fields are little-endian, unless otherwise noted

        buf = struct.pack('<?II', numpy.all(data['DataOK']), len(data) - 1, len(data.values()[0]))
        for i in range(len(data)):
            if data.keys()[i] == 'DataOK':
                continue
            buf += struct.pack('16s', data.keys()[i].encode())
            assert(len(data.values()[i]) == len(data.values()[0]))
            buf += struct.pack('<%di' % len(data.values()[i]), *data.values()[i])

        buf = struct.pack('!q', len(buf)) + buf

        return buf

