#!/usr/bin/env python

from spt3g import core, dfmux, gcp
import socket, argparse, os

parser = argparse.ArgumentParser(description='''
Record dfmux data to disk.

Two modes of operation are suppored.

The standard mode is to supply a hardware map from which to determine which
IceBoards to communicate with.

An alternative mode is to also supply a list of IceBoard serial numbers.  In
that case, the HWM is not loaded from disk (saving as much as several minutes,
depending on the size of the HWM), and the serial numbers are used directly to
initialize the DAQ.  The hardware map file path is instead used to find a
nominal_online_cal.g3 file in the same directory as the YAML file, which is
injected into the DAQ stream for teeing to live visualization tools like
Lyrebird.

Optional auxiliary modules such as the SPT calibrator or the Polarbear WHWP
encoder stream are not enabled by default.

Two methods of error monitoring are also available: enabling the syslog logger
sends all output log messages to a central rsyslog database (e.g. as set up by
SPT), and enabling the GCP watchdog module pings the SPT pager server
periodically to indicate normal operation.
'''
)
parser.add_argument('hardware_map', metavar='/path/to/hwm.yaml', help='Path to hardware map YAML file')
parser.add_argument('boards', nargs='*', metavar='serial', help='IceBoard serial(s) from which to collect data')
parser.add_argument('output', metavar='/path/to/files/', help='Directory in which to place output files')

parser.add_argument('-v', '--verbose', action='store_true', help='Verbose mode (print all frames)')
parser.add_argument('-u', '--udp', action='store_true', help='Use UDP multicast instead of SCTP')
parser.add_argument('--daq-threads', default=None, type=int, help='Number of listener threads to use (applies only for SCTP)')
parser.add_argument('--max-file-size', default=1024, type=int, help='Maximum file size in MB (default 1024)')
parser.add_argument('--calibrator', action='store_true', help='Enable SPT calibrator DAQ')
parser.add_argument('--whwp', action='store_true', help='Enable Polarbear HWP DAQ')
parser.add_argument('--syslog', action='store_true', help='Enable syslog logger')
parser.add_argument('--watchdog', action='store_true', help='Enable GCP watchdog ping')
args = parser.parse_args()

# Tee log messages to both log file and GCP socket
console_logger = core.G3PrintfLogger()
console_logger.timestamps = True # Make sure to get timestamps in the logs

if args.syslog:
    import syslog
    core.G3Logger.global_logger = core.G3MultiLogger(
        [console_logger, core.G3SyslogLogger("dfmuxdaq: ", syslog.LOG_USER)]
    )
else:
    core.G3Logger.global_logger = console_logger

core.set_log_level(core.G3LogLevel.LOG_ERROR, 'DfMuxBuilder')
core.set_log_level(core.G3LogLevel.LOG_ERROR, 'DfMuxCollector')

args.hardware_map = os.path.realpath(args.hardware_map)

if not len(args.boards):
    # If the input is a hardware map path, import the HWM and
    # extract the list of boards from it

    core.log_notice('Initializing hardware map and boards',
                    unit='Data Acquisition')
    import pydfmux
    hwm = pydfmux.load_session(open(args.hardware_map, 'r'))['hardware_map']
    boards = hwm.query(pydfmux.Dfmux)
    boards.resolve()
    boards = boards.serial

else:
    # Otherwise assume the input is a list of board serials
    core.log_notice('Acquiring hardware map information from boards',
                    unit='Data Acquisition')
    hwm = None
    boards = ['%04d' % (int(b)) for b in args.boards]

core.log_notice('Beginning data acquisition', unit='Data Acquisition')
# Set up DfMux consumer
pipe = core.G3Pipeline()
builder = dfmux.DfMuxBuilder([int(board) for board in boards])

if args.udp:
    # Set up listeners per thread and point them at the event builder
    # Get the local IP(s) to use to connect to the boards by opening test
    # connections. Using a set rather than a list deduplicates the results.
    local_ips = {}
    for board in boards:
        testsock = socket.create_connection(('iceboard' + board + '.local', 80))
        local_ip = testsock.getsockname()[0]
        if local_ip not in local_ips:
            local_ips[local_ip] = set()
        local_ips[local_ip].add(int(board))
        testsock.close()
    core.log_notice('Creating listeners for %d boards on interfaces: %s' %
                    (len(boards), ', '.join(local_ips.keys())),
                    unit='Data Acquisition')

    # Set up listeners per network segment and point them at the event builder
    collectors = [dfmux.DfMuxCollector(ip, builder, list(local_boards))
                  for ip, local_boards in local_ips.items()]
else:
    iceboard_hosts = ['iceboard' + board + '.local' for board in boards]
    if args.daq_threads is None:
        args.daq_threads = len(iceboard_hosts)
    if args.daq_threads > len(iceboard_hosts):
        args.daq_threads = len(iceboard_hosts)
    collectors = [dfmux.DfMuxCollector(builder, iceboard_hosts[i::args.daq_threads])
                  for i in range(args.daq_threads)]
pipe.Add(builder)

# Catch errors if samples have become misaligned. Don't even bother processing
# data involving a significant (N-2) reduction from the boards we should have.
nboards = len(boards)
n_badpkts = 0
n_goodpkts = 0
logtweaked = False
def yellanddrop(fr):
    global n_badpkts
    global n_goodpkts
    global logtweaked

    if fr.type != core.G3FrameType.Timepoint:
        return

    if len(fr['DfMux']) < nboards-2:
        if n_badpkts > 0 and n_badpkts % 100 == 0:
            core.log_error('Only %d/%d boards (%s) responding for %d samples -- check for sample misalignment. Temporarily suppressing DfMuxBuilder logging and disabling data archiving.' % (len(fr['DfMux']), nboards, ', '.join([str(k) for k in fr['DfMux'].keys()]), n_badpkts), unit='Data Acquisition')
            # Turn up the threshold on DfMuxBuilder to prevent flooding the console
            core.set_log_level(core.G3LogLevel.LOG_ERROR, 'DfMuxBuilder')
            logtweaked = True
        n_badpkts += 1
        n_goodpkts = 0
        return []
    else:
        n_goodpkts += 1
        if n_goodpkts > 5 and logtweaked:
            # Turn the threshold back down
            core.set_log_level(core.G3LogLevel.LOG_NOTICE, 'DfMuxBuilder')
            core.log_notice('Gross board misalignment resolved. Re-enabling DfMuxBuilder logging and data archiving.', unit='Data Acquisition')
            logtweaked = False
            n_badpkts = 0
pipe.Add(yellanddrop)

# Collect housekeeping when GCP asks for it and send it back to GCP when desired
pipe.Add(gcp.GCPSignalledHousekeeping)
pipe.Add(dfmux.HousekeepingConsumer)
pipe.Add(gcp.GCPHousekeepingTee)

def WaitForWiring(frame):
    if frame.type == core.G3FrameType.Wiring:
        core.set_log_level(core.G3LogLevel.LOG_NOTICE, 'DfMuxBuilder')
        core.set_log_level(core.G3LogLevel.LOG_NOTICE, 'DfMuxCollector')
        core.log_notice('Got a wiring frame, ready for data acquisition.', unit='Data Acquisition')
pipe.Add(WaitForWiring)

# For visualization, add nominal pointing
if hwm is None:
    # Get the bolometer properties map from disk (written separately by pydfmux)
    # whenever a new wiring frame is received.
    def BolometerPropertiesInjector(frame):
        if frame.type == core.G3FrameType.Wiring:
            nchan = len(frame['WiringMap'].keys())
            core.log_notice("Collecting data from %d mapped channels." % (nchan),
                            unit='Data Acquisition')
            bpm = os.path.join(os.path.dirname(args.hardware_map), 'nominal_online_cal.g3')
            try:
                if not os.path.exists(bpm):
                    raise IOError('Missing file %s' % bpm)
                fr = list(core.G3File(bpm))[0]
                return [frame, fr]
            except Exception as e:
                core.log_warn('Error loading BolometerPropertiesMap: %s' % (str(e)),
                              unit='Data Acquisition')
    pipe.Add(BolometerPropertiesInjector)
else:
    # Otherwise inject using the loaded HWM
    def CheckWiring(frame):
        if frame.type == core.G3FrameType.Wiring:
            nchan = len(frame['WiringMap'].keys())
            core.log_notice("Collecting data from %d mapped channels." % (nchan),
                            unit='Data Acquisition')
    pipe.Add(CheckWiring)
    pipe.Add(dfmux.PyDfMuxBolometerPropertiesInjector, pydfmux_hwm=hwm)

# Provide a tee for realtime visualization before possible buffering
# in the calibrator DAQ code
pipe.Add(core.G3ThrottledNetworkSender, hostname='*', port=5451, max_queue_size=1000)

# Collect data from the SPT chopped calibration source
if args.calibrator:
    from spt3g import auxdaq
    pipe.Add(auxdaq.CalibratorDAQ)

# Collect data from the Polarbear warm half wave plate angular encoder
if args.whwp:
    from spt3g import whwp
    pipe.Add(whwp.WHWPConsumer)

# Issue a periodic watchdog ping to the SPT pager system
if args.watchdog:
    pipe.Add(gcp.DAQWatchdog, calibrator=args.calibrator)

if args.verbose:
    pipe.Add(core.Dump)

# Provide a tee for other software to collect frames
pipe.Add(core.G3NetworkSender, hostname='*', port=5353, max_queue_size=1000)

def filename(frame, seq):
    if 'EventHeader' in frame:
        t = frame['EventHeader']
    else:
        t = core.G3Time.Now()
    return os.path.join(args.output, t.GetFileFormatString() + '.g3')
pipe.Add(core.G3MultiFileWriter, filename=filename, size_limit=args.max_file_size*1024*1024)

for collector in collectors:
    collector.Start()
pipe.Run(profile=True)

# Shut everything down
for collector in collectors:
    collector.Stop()

# C++ global destructor runs after python has closed, sometimes
# Setting the logger to None like this explicitly avoids segfaults
core.G3Logger.global_logger = None
