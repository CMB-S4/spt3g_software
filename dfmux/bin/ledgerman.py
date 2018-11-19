#!/usr/bin/env python

from spt3g import core, dfmux
import socket, argparse
import numpy as np

parser = argparse.ArgumentParser(description='Record dfmux data to a NetCDF file', prog='ledgerman')
parser.add_argument('hardware_map', metavar='/path/to/hwm.yaml', help='Path to hardware map YAML file')
parser.add_argument('output', metavar='output.nc', help='Path to output NetCDF file')

parser.add_argument('-v', dest='verbose', action='store_true',
                    help='Verbose mode (print all frames)')
parser.add_argument('-a', dest='align', action='store_true',
                    help='Align sampling. This has to happen once, but'
                    ' will break any existing DAN loops when run.')
parser.add_argument('-s', dest='system_time', action='store_true',
                    help='Replace board time with system time when data'
                    ' received. Useful if your board is in IRIG_TEST mode.')
parser.add_argument('-P', dest='physnames', action='store_true',
                    help='Use physical bolometer names rather than channel names')
parser.add_argument('-p', dest='pathstring', action='store', default=[], nargs='+',
                    help='Only record channels that match channel path string,'
                    ' e.g. 005/5/2/3/* saves data for all bolometers on crate'
                    ' 005, slot 5, mezzanine 2, module 3.')
parser.add_argument('-b', dest='state', action='store', default=[], nargs='+',
                    help='Only record bolometers that have a state that'
                    ' matches the supplied string(s), e.g. overbiased tuned')
args = parser.parse_args()


# Import pydfmux later since it can take a while
import pydfmux

core.log_notice('Initializing hardware map and boards', unit='Ledgerman')
hwm = pydfmux.load_session(open(args.hardware_map, 'r'))['hardware_map']
if hwm.query(pydfmux.IceCrate).count() > 0:
    hwm.query(pydfmux.IceCrate).resolve()

# make sure that the hardware map is consistent with what's on the IceBoards
if args.state:
    hwm.query(pydfmux.Bolometer).load_bolo_states()

if args.align:
    core.log_notice('Aligning board sampling, this will break any existing DAN loops!',
                    unit='Ledgerman')
    hwm.query(pydfmux.IceBoard).set_fir_stage(6)
    hwm.query(pydfmux.IceBoard).align_sampling()

core.log_notice('Beginning data acquisition', unit='Ledgerman')

# get board serial numbers only for the channels that we are recording
if args.pathstring:
    chan_map_query = hwm.channel_maps_from_pstring(args.pathstring)
else:
    chan_map_query = hwm.query(pydfmux.ChannelMapping)
if args.state:
    chan_map_query = chan_map_query.join(pydfmux.ChannelMapping, pydfmux.Bolometer).filter(pydfmux.Bolometer.state._in(args.state))
serial_list = np.unique(np.array([cm.iceboard.serial for cm in chan_map_query]))

# Set up DfMux consumer
pipe = core.G3Pipeline()
builder = dfmux.DfMuxBuilder([int(serial) for serial in serial_list])

# Get the local IP(s) to use to connect to the boards by opening test
# connections. Using a set rather than a list deduplicates the results.
local_ips = set()
for board in hwm.query(pydfmux.core.dfmux.IceBoard):
    testsock = socket.create_connection(('iceboard' + board.serial + '.local', 80))
    local_ips.add(testsock.getsockname()[0])
    testsock.close()
core.log_notice('Creating listeners for %d boards on interfaces: %s' % (hwm.query(pydfmux.core.dfmux.IceBoard).count(), ', '.join(local_ips)), unit='Ledgerman')

# Build mapping dictionary for old (64x) firmware
v2_mapping = {'iceboard' + str(serial) + '.local': int(serial) for serial in serial_list}

# Set up listeners per network segment and point them at the event builder
collectors = [dfmux.DfMuxCollector(ip, builder, v2_mapping) for ip in local_ips]
pipe.Add(builder)

# Insert current hardware map into data stream. This is critical to get the
# board ID -> IP mapping needed to do anything useful with the data
pipe.Add(dfmux.PyDfMuxHardwareMapInjector, pydfmux_hwm=hwm, pathstring=args.pathstring, state=args.state)

if args.physnames:
    from spt3g import calibration
    # Swap the normal bolometer names for their physical_name equivalents,
    # which can make laboratory tasks simpler -- we don't care about
    # long-term data archiving here.
    pipe.Add(dfmux.PyDfMuxBolometerPropertiesInjector, pydfmux_hwm=hwm)

    class SwapPhysLogicalNames(object):
        def __init__(self):
            self.calframe = None
            self.wiringframe = None
            self.sent = False

        def __call__(self, frame):
            if self.sent:
                return
            if frame.type == core.G3FrameType.Wiring:
                self.wiringframe = frame
                return []
            if frame.type == core.G3FrameType.Calibration:
                self.calframe = frame
                return []
            if self.calframe is None or self.wiringframe is None:
                raise Exception('Data before wiring and cal frames!')
            w = self.wiringframe['WiringMap']
            del self.wiringframe['WiringMap']
            c = self.calframe['NominalBolometerProperties']

            new_wiring = dfmux.DfMuxWiringMap()
            for k, d in w.iteritems():
                new_wiring[c[k].physical_name] = d
            self.wiringframe['WiringMap'] = new_wiring

            self.sent = True
            return [self.wiringframe, frame]
    pipe.Add(SwapPhysLogicalNames)

if args.system_time:
    def sub_system_time(frame):
        if frame.type != core.G3FrameType.Timepoint:
            return
        del frame['EventHeader']
        frame['EventHeader'] = core.G3Time.Now()
    pipe.Add(sub_system_time)

if args.verbose:
    pipe.Add(core.Dump)


pipe.Add(dfmux.NetCDFDump, filename=args.output)

for collector in collectors:
    collector.Start()
pipe.Run()

