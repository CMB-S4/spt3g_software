#!/usr/bin/env python

from spt3g import core, dfmux, gcp, auxdaq
import socket, argparse, os

parser = argparse.ArgumentParser(description='Record dfmux data to disk')
parser.add_argument('hardware_map', metavar='/path/to/hwm.yaml', help='Path to hardware map YAML file')
parser.add_argument('output', metavar='/path/to/files/', help='Directory in which to place output files')

parser.add_argument('-v', dest='verbose', action='store_true', help='Verbose mode (print all frames)')
parser.add_argument('--max_file_size', dest='max_file_size', default=1024, help='Maximum file size in MB (default 1024)')
args = parser.parse_args()

# Import pydfmux later since it can take a while
import pydfmux

print('Initializing hardware map and boards')
hwm = pydfmux.load_session(open(args.hardware_map, 'r'))['hardware_map']
if hwm.query(pydfmux.core.dfmux.IceCrate).count() > 0:
        hwm.query(pydfmux.core.dfmux.IceCrate).resolve()

print('Beginning data acquisition')
# Set up DfMux consumer
pipe = core.G3Pipeline()
builder = dfmux.DfMuxBuilder([int(board.serial) for board in hwm.query(pydfmux.IceBoard)])

# Get the local IP(s) to use to connect to the boards by opening test
# connections. Using a set rather than a list deduplicates the results.
local_ips = set()
for board in hwm.query(pydfmux.core.dfmux.IceBoard):
    testsock = socket.create_connection(('iceboard' + board.serial + '.local', 80))
    local_ips.add(testsock.getsockname()[0])
    testsock.close()
print('Creating listeners for %d boards on interfaces: %s' % (hwm.query(pydfmux.core.dfmux.IceBoard).count(), ', '.join(local_ips)))

# Set up listeners per network segment and point them at the event builder
collectors = [dfmux.DfMuxCollector(ip, builder) for ip in local_ips]
pipe.Add(builder)

# Catch errors if samples have become misaligned. Don't even bother processing
# data involving a significant (N-2) reduction from the boards we should have.
nboards = hwm.query(pydfmux.IceBoard).count()
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

# Insert current hardware map into data stream. This is critical to get the
# board ID -> IP mapping needed to do anything useful with the data
pipe.Add(dfmux.PyDfMuxHardwareMapInjector, pydfmux_hwm=hwm)

# For visualization, add nominal pointing
pipe.Add(dfmux.PyDfMuxBolometerPropertiesInjector, pydfmux_hwm=hwm)

# Collect housekeeping when GCP asks for it and send it back to GCP when desired
pipe.Add(gcp.GCPSignalledHousekeeping)
pipe.Add(dfmux.HousekeepingConsumer)
pipe.Add(gcp.GCPHousekeepingTee)

# Provide a tee for realtime visualization before possible buffering
# in the calibrator DAQ code
pipe.Add(core.G3ThrottledNetworkSender, hostname='*', port=5451, max_queue_size=1000)

# Collect data from the chopped calibration source
pipe.Add(auxdaq.CalibratorDAQ)

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
pipe.Run()

